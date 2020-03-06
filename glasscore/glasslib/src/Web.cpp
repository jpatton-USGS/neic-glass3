#include "Web.h"
#include <json.h>
#include <logger.h>
#include <geo.h>
#include <cmath>
#include <algorithm>
#include <string>
#include <tuple>
#include <memory>
#include <utility>
#include <vector>
#include <fstream>
#include <thread>
#include <mutex>
#include <ctime>
#include <limits>
#include <map>
#include <sstream>
#include "Glass.h"
#include "Pick.h"
#include "Node.h"
#include "SiteList.h"
#include "Site.h"

namespace glasscore {

// constants
constexpr double CWeb::k_dAzimuthTaperDefault;
constexpr double CWeb::k_dDepthResolutionUndefined;
constexpr double CWeb::k_dFibonacciRatio;
const int CWeb::k_iNodeLatitudeIndex;
const int CWeb::k_iNodeLongitudeIndex;
const int CWeb::k_iNodeDepthIndex;
constexpr double CWeb::k_dMinimumMaxNodeDepth;

// site sorting function
// Compares nodal distance for nearest site assignment
bool sortSite(const std::pair<double, std::shared_ptr<CSite>> &lhs,
				const std::pair<double, std::shared_ptr<CSite>> &rhs) {
	// compare
	if (lhs.first < rhs.first) {
		return (true);
	}

	// lhs > rhs
	return (false);
}

// ---------------------------------------------------------CWeb
CWeb::CWeb(int numThreads, int sleepTime, int checkInterval)
		: glass3::util::ThreadBaseClass("Web", sleepTime, numThreads,
										checkInterval) {
	clear();

	// start up the threads
	start();
}

// ---------------------------------------------------------CWeb
CWeb::CWeb(std::string name, double thresh, int numDetect, int numNucleate,
			int resolution, bool update, bool save,
			std::shared_ptr<traveltime::CTravelTime> firstTrav,
			std::shared_ptr<traveltime::CTravelTime> secondTrav, int numThreads,
			int sleepTime, int checkInterval, double aziTap, double maxDep)
		: glass3::util::ThreadBaseClass("Web", sleepTime, numThreads,
										checkInterval) {
	clear();

	initialize(name, thresh, numDetect, numNucleate, resolution, update, save,
				firstTrav, secondTrav, aziTap, maxDep);

	// start up the threads
	start();
}

// ---------------------------------------------------------~CWeb
CWeb::~CWeb() {
}

// ---------------------------------------------------------clear
void CWeb::clear() {
	std::lock_guard<std::recursive_mutex> webGuard(m_WebMutex);
	m_iNumStationsPerNode = 10;
	m_iNucleationDataCountThreshold = 5;
	m_dNucleationStackThreshold = 2.5;
	m_dNodeResolution = 100;
	m_dDepthResolution = k_dDepthResolutionUndefined;
	m_sName = "UNDEFINED";
	m_pSiteList = NULL;
	m_bUpdate = false;
	m_bSaveGrid = false;
	m_dAzimuthTaper = k_dAzimuthTaperDefault;
	m_dMaxDepth = CGlass::k_dMaximumDepth;

	// clear out all the nodes in the web
	try {
		m_vNodeMutex.lock();
		for (auto &node : m_vNode) {
			node->clear();
		}
		m_vNode.clear();
	} catch (...) {
		// ensure the vNode mutex is unlocked
		m_vNodeMutex.unlock();

		throw;
	}
	m_vNodeMutex.unlock();

	// clear the network filter
	m_vNetworksFilter.clear();
	m_vSitesFilter.clear();
	m_bUseOnlyTeleseismicStations = false;

	// clear sites
	try {
		m_vSiteMutex.lock();
		vWebSites.clear();
	} catch (...) {
		// ensure the vSite mutex is unlocked
		m_vSiteMutex.unlock();

		throw;
	}
	m_vSiteMutex.unlock();

	m_pNucleationTravelTime1 = NULL;
	m_pNucleationTravelTime2 = NULL;

	// reset zonestats info
	m_pZoneStats = NULL;
	m_sZoneStatsFileName.clear();

	// reset quality filter
	m_dQualityFilter = -1.0;
}

// ---------------------------------------------------------initialize
bool CWeb::initialize(std::string name, double thresh, int numDetect,
						int numNucleate, int resolution, bool update, bool save,
						std::shared_ptr<traveltime::CTravelTime> firstTrav,
						std::shared_ptr<traveltime::CTravelTime> secondTrav,
						double aziTap, double maxDep) {
	std::lock_guard<std::recursive_mutex> webGuard(m_WebMutex);

	m_sName = name;
	m_dNucleationStackThreshold = thresh;
	m_iNumStationsPerNode = numDetect;
	m_iNucleationDataCountThreshold = numNucleate;
	m_dNodeResolution = resolution;
	m_bUpdate = update;
	m_bSaveGrid = save;

	m_travelTimeMutex.lock();
	m_pNucleationTravelTime1 = firstTrav;
	m_pNucleationTravelTime2 = secondTrav;
	m_travelTimeMutex.unlock();

	m_dAzimuthTaper = aziTap;
	m_dMaxDepth = maxDep;

	m_tLastUpdated = -1;

	// done
	return (true);
}

// -------------------------------------------------------receiveExternalMessage
bool CWeb::receiveExternalMessage(std::shared_ptr<json::Object> com) {
	// null check json
	if (com == NULL) {
		glass3::util::Logger::log(
				"error",
				"CWeb::receiveExternalMessage: NULL json communication.");
		return (false);
	}

	// we only care about messages with a string Cmd key.
	if (com->HasKey("Cmd")
			&& ((*com)["Cmd"].GetType() == json::ValueType::StringVal)) {
		// dispatch to appropriate function based on Cmd value
		json::Value v = (*com)["Cmd"].ToString();

		// generate a Global grid
		if (v == "Global") {
			return (generateGlobalGrid(com));
		}

		// generate a regional or local grid
		if (v == "Grid") {
			return (generateLocalGrid(com));
		}

		// generate an explicitly defined grid
		if (v == "GridExplicit") {
			return (generateExplicitGrid(com));
		}
	}

	// this communication was not handled
	return (false);
}

// ---------------------------------------------------------generateGlobalGrid
bool CWeb::generateGlobalGrid(std::shared_ptr<json::Object> gridConfiguration) {
	glass3::util::Logger::log("debug", "CWeb::generateGlobalGrid");

	// nullchecks
	// check json
	if (gridConfiguration == NULL) {
		glass3::util::Logger::log(
				"error", "CWeb::generateGlobalGrid: NULL json configuration.");
		return (false);
	}

	char sLog[glass3::util::Logger::k_nMaxLogEntrySize];
	std::vector<double> depthLayerArray;
	int numDepthLayers = 0;

	// load basic (common) grid configuration
	if (loadGridConfiguration(gridConfiguration) == false) {
		return (false);
	}

	// list of depth layers to generate for this grid
	if (((*gridConfiguration).HasKey("DepthLayers"))
			&& ((*gridConfiguration)["DepthLayers"].GetType()
					== json::ValueType::ArrayVal)) {
		json::Array zarray = (*gridConfiguration)["DepthLayers"].ToArray();
		for (auto v : zarray) {
			if (v.GetType() == json::ValueType::DoubleVal) {
				double depth = v.ToDouble();
				depthLayerArray.push_back(depth);
			}
		}
		numDepthLayers = static_cast<int>(depthLayerArray.size());
	} else {
		glass3::util::Logger::log(
				"error",
				"CWeb::generateGlobalGrid: Missing required DepthLayers Array.");
		return (false);
	}

	// calculate the number of nodes from the desired resolution
	// using a function that was EMPIRICALLY determined via using different
	// numNode values and computing the average resolution from a node to
	// the nearest other 6 nodes. The spreadsheet used to calculate this function
	// is located in ** NodesToResoultionCalculations.xlsx **
	// The intention is to calculate the number of nodes to ensure the desired
	// node resolution when the grid is generated below
	int numNodes = 5.3E8 * std::pow(getNodeResolution(), -1.965);

	// should have an odd number of nodes (see paper named below)
	if ((numNodes % 2) == 0) {
		numNodes += 1;
	}

	snprintf(sLog, sizeof(sLog),
				"CWeb::generateGlobalGrid: Calculated numNodes:%d;", numNodes);
	glass3::util::Logger::log(sLog);

	// lock the site list while getting site list
	while ((m_vSiteMutex.try_lock() == false) &&
					(getTerminate() == false)) {
		// update thread status
		setThreadHealth(true);

		// wait a little while
		std::this_thread::sleep_for(
				std::chrono::milliseconds(getSleepTime()));
	}

	// grab a copy of the site vector for our use
	std::vector<std::pair<double, std::shared_ptr<CSite>>>
		sites(vWebSites);

	m_vSiteMutex.unlock();

	// create / open gridfile for saving
	std::ofstream outfile;
	std::ofstream outstafile;
	if (getSaveGrid()) {
		std::string filename = m_sName + "_gridfile.txt";
		outfile.open(filename, std::ios::out);
		outfile << "generateLocalGrid,NodeID,NodeLat,NodeLon,NodeDepth" << "\n";

		filename = m_sName + "_gridstafile.txt";
		outstafile.open(filename, std::ios::out);
		outstafile << "NodeID,[StationSCNL;StationLat;StationLon;StationRad],"
					<< "\n";
	}

	// Generate equally spaced grid of nodes over the globe (more or less)
	// Follows Paper (Gonzalez, 2010) Measurement of Areas on a Sphere Using
	// Fibonacci and Latitude Longitude Lattices
	// std::vector<std::pair<double, double>> vVert;
	int iNodeCount = 0;
	int numSamples = (numNodes - 1) / 2;

	for (int i = (-1 * numSamples); i <= numSamples; i++) {
		double aLat = std::asin((2 * i) / ((2.0 * numSamples) + 1))
				* (180.0 / glass3::util::GlassMath::k_Pi);
		double aLon = fmod(i, k_dFibonacciRatio) * (360.0 / k_dFibonacciRatio);

		// longitude bounds check
		if (aLon < glass3::util::Geo::k_MinimumLongitude) {
			aLon += glass3::util::Geo::k_LongitudeWrap;
		}
		if (aLon > glass3::util::Geo::k_MaximumLongitude) {
			aLon -= glass3::util::Geo::k_LongitudeWrap;
		}

		// use zonestats to get the max depth for this node, if we have
		// zonestats available, otherwise default to the configured
		// max depth for the grid
		double dMaxNodeDepth = m_dMaxDepth;
		// bool bReachedMaxDepth = false;
		if (m_pZoneStats != NULL) {
			double aDepth = m_pZoneStats->getMaxDepthForLatLon(aLat, aLon);
			if (aDepth != m_pZoneStats->depthInvalid) {
				dMaxNodeDepth = aDepth;
			}
		}

		// for each depth
		for (auto z : depthLayerArray) {
			// check to see if Z is below the maximum depth
			// if (z >= std::max(dMaxNodeDepth, k_dMinimumMaxNodeDepth)) {
			// bReachedMaxDepth = true;
			// z = dMaxNodeDepth;
			// }

			if (z > std::max(dMaxNodeDepth, k_dMinimumMaxNodeDepth)) {
				break;
			}

			// sort site list for this vertex
			sites = sortSiteListForNode(aLat, aLon, z, sites);

			// it would make a certain amount of sense here, to track
			// the depth delta between this node and the vertically
			// adjacent ones (ones above and below it), and save
			// the max of those as the Vertical Resolution, which
			// could be used to set boundaries during nucleation,
			// and POSSIBLY also constrain the solution post-nucleation,
			// during initial location.
			// create node
			std::shared_ptr<CNode> node = generateNode(aLat, aLon, z,
														getNodeResolution(), sites);

			// if we got a valid node, add it
			if (addNode(node) == true) {
				iNodeCount++;

				// write node to generateLocalGrid file
				if (getSaveGrid()) {
					outfile << m_sName << "," << node->getID() << ","
							<< std::to_string(aLat) << ","
							<< std::to_string(aLon) << "," << std::to_string(z)
							<< "\n";

					// write to station file
					outstafile << node->getSitesString();
				}
			}  // end if addNode()

			// if (bReachedMaxDepth) {
			// break;
			// }
		}  // end for each depth in depthLayerArray
	}  // end for each sample

	// close generateLocalGrid file
	if (getSaveGrid()) {
		outfile.close();
		outstafile.close();
	}

	std::string phases = "";
	if (m_pNucleationTravelTime1 != NULL) {
		std::lock_guard<std::recursive_mutex> ttGuard(m_travelTimeMutex);
		phases += m_pNucleationTravelTime1->m_sPhase;
	}
	if (m_pNucleationTravelTime2 != NULL) {
		std::lock_guard<std::recursive_mutex> ttGuard(m_travelTimeMutex);
		phases += ", " + m_pNucleationTravelTime2->m_sPhase;
	}

	// log grid info
	snprintf(
			sLog, sizeof(sLog),
			"CWeb::generateGlobalGrid sName:%s Phase(s):%s; nZ:%d; resol:%.2f;"
			" nDetect:%d; nNucleate:%d; dThresh:%.2f; vNetFilter:%d;"
			" vSitesFilter:%d; bUseOnlyTeleseismicStations:%d; iNodeCount:%d;",
			m_sName.c_str(), phases.c_str(), numDepthLayers,
			getNodeResolution(), getNumStationsPerNode(),
			getNucleationDataCountThreshold(), getNucleationStackThreshold(),
			static_cast<int>(m_vNetworksFilter.size()),
			static_cast<int>(m_vSitesFilter.size()),
			static_cast<int>(m_bUseOnlyTeleseismicStations), iNodeCount);
	glass3::util::Logger::log("info", sLog);

	// success
	return (true);
}

// ---------------------------------------------------------generateLocalGrid
bool CWeb::generateLocalGrid(std::shared_ptr<json::Object> gridConfiguration) {
	glass3::util::Logger::log("debug", "CWeb::generateLocalGrid");

	// nullchecks
	// check json
	if (gridConfiguration == NULL) {
		glass3::util::Logger::log(
				"error", "CWeb::generateLocalGrid: NULL json configuration.");
		return (false);
	}

	char sLog[glass3::util::Logger::k_nMaxLogEntrySize];
	double lat = 0;
	double lon = 0;
	int rows = 0;
	int cols = 0;
	std::vector<double> depthLayerArray;
	int numDepthLayers = 0;

	// load basic (common) grid configuration
	if (loadGridConfiguration(gridConfiguration) == false) {
		return (false);
	}

	// Latitude of the center point of this generateLocalGrid
	if (((*gridConfiguration).HasKey("CenterLatitude"))
			&& ((*gridConfiguration)["CenterLatitude"].GetType()
					== json::ValueType::DoubleVal)) {
		lat = (*gridConfiguration)["CenterLatitude"].ToDouble();
	} else {
		glass3::util::Logger::log(
				"error",
				"CWeb::generateLocalGrid: Missing required CenterLatitude Key.");
		return (false);
	}

	// Longitude of the center point of this generateLocalGrid
	if (((*gridConfiguration).HasKey("CenterLongitude"))
			&& ((*gridConfiguration)["CenterLongitude"].GetType()
					== json::ValueType::DoubleVal)) {
		lon = (*gridConfiguration)["CenterLongitude"].ToDouble();
	} else {
		glass3::util::Logger::log(
				"error", "CWeb::generateLocalGrid: Missing required Lon Key.");
		return (false);
	}

	// list of depth layers to generate for this generateLocalGrid
	if (((*gridConfiguration).HasKey("DepthLayers"))
			&& ((*gridConfiguration)["DepthLayers"].GetType()
					== json::ValueType::ArrayVal)) {
		json::Array jsonLayers = (*gridConfiguration)["DepthLayers"].ToArray();
		for (auto aLayer : jsonLayers) {
			if (aLayer.GetType() == json::ValueType::DoubleVal) {
				depthLayerArray.push_back(aLayer.ToDouble());
			}
		}
		numDepthLayers = static_cast<int>(depthLayerArray.size());
	} else {
		glass3::util::Logger::log(
				"error",
				"CWeb::generateLocalGrid: Missing required DepthLayers Array.");
		return (false);
	}

	// the number of rows in this generateLocalGrid
	if (((*gridConfiguration).HasKey("NumberOfRows"))
			&& ((*gridConfiguration)["NumberOfRows"].GetType()
					== json::ValueType::IntVal)) {
		rows = (*gridConfiguration)["NumberOfRows"].ToInt();
	} else {
		glass3::util::Logger::log(
				"error",
				"CWeb::generateLocalGrid: Missing required NumberOfRows Key.");
		return (false);
	}

	// the number of columns in this generateLocalGrid
	if (((*gridConfiguration).HasKey("NumberOfColumns"))
			&& ((*gridConfiguration)["NumberOfColumns"].GetType()
					== json::ValueType::IntVal)) {
		cols = (*gridConfiguration)["NumberOfColumns"].ToInt();
	} else {
		glass3::util::Logger::log(
				"error",
				"CWeb::generateLocalGrid: Missing required NumberOfColumns Key.");
		return (false);
	}

	// Generate initial node values
	// compute latitude distance in geographic degrees by converting
	// the provided resolution in kilometers to degrees
	// NOTE: Hard coded conversion factor
	double latDistance = getNodeResolution() / glass3::util::Geo::k_DegreesToKm;

	// compute the longitude distance in geographic degrees by
	// dividing the latitude distance by the cosine of the center latitude
	double lonDistance = latDistance
			/ cos(glass3::util::GlassMath::k_DegreesToRadians * lat);

	// compute the middle row index
	int irow0 = rows / 2;

	// compute the middle column index
	int icol0 = cols / 2;

	// compute the maximum latitude using the provided center latitude,
	// middle row index, and the latitude distance
	double lat0 = lat + (irow0 * latDistance);

	// compute the minimum longitude using the provided center longitude,
	// middle column index, and the longitude distance
	double lon0 = lon - (icol0 * lonDistance);

	// lock the site list while getting site list
	while ((m_vSiteMutex.try_lock() == false) &&
					(getTerminate() == false)) {
		// update thread status
		setThreadHealth(true);

		// wait a little while
		std::this_thread::sleep_for(
				std::chrono::milliseconds(getSleepTime()));
	}

	// grab a copy of the site vector for our use
	std::vector<std::pair<double, std::shared_ptr<CSite>>>
		sites(vWebSites);

	m_vSiteMutex.unlock();

	// create / open gridfile for saving
	std::ofstream outfile;
	std::ofstream outstafile;
	if (getSaveGrid()) {
		std::string filename = m_sName + "_gridfile.txt";
		outfile.open(filename, std::ios::out);
		outfile << "generateLocalGrid,NodeID,NodeLat,NodeLon,NodeDepth" << "\n";

		filename = m_sName + "_gridstafile.txt";
		outstafile.open(filename, std::ios::out);
		outstafile << "NodeID,[StationSCNL;StationLat;StationLon;StationRad],"
					<< "\n";
	}

	// init node count
	int iNodeCount = 0;

	// generate grid
	// for each row
	for (int irow = 0; irow < rows; irow++) {
		// compute the current row latitude by subtracting
		// the row index multiplied by the latitude distance from the
		// maximum latitude
		double latrow = lat0 - (irow * latDistance);

		// log the latitude
		snprintf(sLog, sizeof(sLog), "LatRow:%.2f", latrow);
		glass3::util::Logger::log("debug", sLog);
		// for each col
		for (int icol = 0; icol < cols; icol++) {
			// compute the current column longitude by adding the
			// column index multiplied by the longitude distance to the
			// minimum longitude
			double loncol = lon0 + (icol * lonDistance);

			// use zonestats to get the max depth for this node, if we have
			// zonestats available, otherwise default to the configured
			// max depth for the grid
			double dMaxNodeDepth = m_dMaxDepth;
			if (m_pZoneStats != NULL) {
				double aDepth = m_pZoneStats->getMaxDepthForLatLon(latrow,
																	loncol);
				if (aDepth != m_pZoneStats->depthInvalid) {
					dMaxNodeDepth = aDepth;
				}
			}

			// for each depth at this generateLocalGrid point
			for (auto z : depthLayerArray) {
				// check to see if Z is below the maximum depth
				if (z > std::max(dMaxNodeDepth, k_dMinimumMaxNodeDepth)) {
					break;
				}

				// sort site list for this generateLocalGrid point
				sites = sortSiteListForNode(latrow, loncol, z, sites);

				// generate this node
				std::shared_ptr<CNode> node = generateNode(
						latrow, loncol, z, getNodeResolution(), sites);

				// if we got a valid node, add it
				if (addNode(node) == true) {
					iNodeCount++;
				}  // end if addNode()

				// write node to generateLocalGrid file
				if (getSaveGrid()) {
					outfile << m_sName << "," << node->getID() << ","
							<< std::to_string(latrow) << ","
							<< std::to_string(loncol) << ","
							<< std::to_string(z) << "\n";

					// write to station file
					outstafile << node->getSitesString();
				}  // end if getSaveGrid()
			}  // end for each depth layer
		}  // end for each lon-column in grid
	}  // end for each lat-row in grid

	// close generateLocalGrid file
	if (getSaveGrid()) {
		outfile.close();
		outstafile.close();
	}

	std::string phases = "";
	if (m_pNucleationTravelTime1 != NULL) {
		std::lock_guard<std::recursive_mutex> ttGuard(m_travelTimeMutex);
		phases += m_pNucleationTravelTime1->m_sPhase;
	}
	if (m_pNucleationTravelTime2 != NULL) {
		std::lock_guard<std::recursive_mutex> ttGuard(m_travelTimeMutex);
		phases += ", " + m_pNucleationTravelTime2->m_sPhase;
	}

	// log local grid info
	snprintf(
			sLog,
			sizeof(sLog),
			"CWeb::generateLocalGrid sName:%s Phase(s):%s; Ranges:Lat(%.2f,%.2f),"
			"Lon:(%.2f,%.2f); nRow:%d; nCol:%d; nZ:%d; resol:%.2f;"
			" nDetect:%d; nNucleate:%d; dThresh:%.2f; vNetFilter:%d;"
			" vSitesFilter:%d;  bUseOnlyTeleseismicStations:%d;"
			" iNodeCount:%d;",
			m_sName.c_str(), phases.c_str(), lat0,
			lat0 - (rows - 1) * latDistance, lon0,
			lon0 + (cols - 1) * lonDistance, rows, cols, numDepthLayers,
			getNodeResolution(), getNumStationsPerNode(),
			getNucleationDataCountThreshold(), getNucleationStackThreshold(),
			static_cast<int>(m_vNetworksFilter.size()),
			static_cast<int>(m_vSitesFilter.size()),
			static_cast<int>(m_bUseOnlyTeleseismicStations), iNodeCount);
	glass3::util::Logger::log("info", sLog);

	// success
	return (true);
}

// ---------------------------------------------------------generateExplicitGrid
bool CWeb::generateExplicitGrid(
		std::shared_ptr<json::Object> gridConfiguration) {
	glass3::util::Logger::log("debug", "CWeb::generateExplicitGrid");

	// nullchecks
	// check json
	if (gridConfiguration == NULL) {
		glass3::util::Logger::log(
				"error",
				"CWeb::generateExplicitGrid: NULL json configuration.");
		return (false);
	}

	char sLog[glass3::util::Logger::k_nMaxLogEntrySize];

	// load basic (common) grid configuration
	if (loadGridConfiguration(gridConfiguration) == false) {
		return (false);
	}

	int nN = 0;
	std::vector<std::vector<double>> nodes;
	double resol = 0;

	// the number of columns in this generateLocalGrid
	if ((gridConfiguration->HasKey("NodeList"))
			&& ((*gridConfiguration)["NodeList"].GetType()
					== json::ValueType::ArrayVal)) {
		json::Array nodeList = (*gridConfiguration)["NodeList"].ToArray();

		for (const auto &val : nodeList) {
			nN++;
		}
		nodes.resize(nN);

		int i = 0;
		for (const auto &val : nodeList) {
			if (val.GetType() != json::ValueType::ObjectVal) {
				glass3::util::Logger::log(
						"error",
						"CWeb::generateExplicitGrid: Bad Node Definition found"
						" in Node List.");
				continue;
			}
			json::Object obj = val.ToObject();
			if (obj["Latitude"].GetType() != json::ValueType::DoubleVal
					|| obj["Longitude"].GetType() != json::ValueType::DoubleVal
					|| obj["Depth"].GetType() != json::ValueType::DoubleVal) {
				glass3::util::Logger::log(
						"error",
						"CWeb::generateExplicitGrid: Node Lat, Lon, or Depth not"
						" a  double in param file.");
				return (false);
			}
			double lat = obj["Latitude"].ToDouble();
			double lon = obj["Longitude"].ToDouble();
			double depth = obj["Depth"].ToDouble();
			nodes[i].resize(3);
			nodes[i][k_iNodeLatitudeIndex] = lat;
			nodes[i][k_iNodeLongitudeIndex] = lon;
			nodes[i][k_iNodeDepthIndex] = depth;
			i++;
		}
	} else {
		glass3::util::Logger::log(
				"error",
				"CWeb::generateExplicitGrid: Missing required NodeList Key.");
		return (false);
	}

	// lock the site list while getting site list
	while ((m_vSiteMutex.try_lock() == false) &&
					(getTerminate() == false)) {
		// update thread status
		setThreadHealth(true);

		// wait a little while
		std::this_thread::sleep_for(
				std::chrono::milliseconds(getSleepTime()));
	}

	// grab a copy of the site vector for our use
	std::vector<std::pair<double, std::shared_ptr<CSite>>>
		sites(vWebSites);

	m_vSiteMutex.unlock();

	// create / open gridfile for saving
	std::ofstream outfile;
	std::ofstream outstafile;
	if (getSaveGrid()) {
		std::string filename = m_sName + "_gridfile.txt";
		outfile.open(filename, std::ios::out);
		outfile << "generateLocalGrid,NodeID,NodeLat,NodeLon,NodeDepth" << "\n";

		filename = m_sName + "_gridstafile.txt";
		outstafile.open(filename, std::ios::out);
		outstafile << "NodeID,StationSCNL,StationLat,StationLon,StationRad"
					<< "\n";
	}

	// init node count
	int iNodeCount = 0;

	// loop through node vector
	for (int i = 0; i < nN; i++) {
		// get lat,lon,depth
		double lat = nodes[i][k_iNodeLatitudeIndex];
		double lon = nodes[i][k_iNodeLongitudeIndex];
		double Z = nodes[i][k_iNodeDepthIndex];

		// sort site list
		sites = sortSiteListForNode(lat, lon, Z, sites);

		// don't do any maxdepth/zonestats checks here, since this grid is
		// explicit

		// create node, note resolution set to 0
		std::shared_ptr<CNode> node = generateNode(lat, lon, Z,
													getNodeResolution(), sites);

		if (addNode(node) == true) {
			iNodeCount++;
		}

		// write to station file
		outstafile << node->getSitesString();
	}

	// close grid file
	if (getSaveGrid()) {
		outfile.close();
		outstafile.close();
	}

	std::string phases = "";
	if (m_pNucleationTravelTime1 != NULL) {
		std::lock_guard<std::recursive_mutex> ttGuard(m_travelTimeMutex);
		phases += m_pNucleationTravelTime1->m_sPhase;
	}
	if (m_pNucleationTravelTime2 != NULL) {
		std::lock_guard<std::recursive_mutex> ttGuard(m_travelTimeMutex);
		phases += ", " + m_pNucleationTravelTime2->m_sPhase;
	}

	// log explicit grid info
	snprintf(
			sLog, sizeof(sLog),
			"CWeb::generateExplicitGrid sName:%s Phase(s):%s; nDetect:%d;"
			" nNucleate:%d; dThresh:%.2f; vNetFilter:%d;"
			" bUseOnlyTeleseismicStations:%d; vSitesFilter:%d; iNodeCount:%d;",
			m_sName.c_str(), phases.c_str(), getNumStationsPerNode(),
			getNucleationDataCountThreshold(), getNucleationStackThreshold(),
			static_cast<int>(m_vNetworksFilter.size()),
			static_cast<int>(m_vSitesFilter.size()),
			static_cast<int>(m_bUseOnlyTeleseismicStations), iNodeCount);
	glass3::util::Logger::log("info", sLog);

	// success
	return (true);
}

// ---------------------------------------------------------loadGridConfiguration
bool CWeb::loadGridConfiguration(
		std::shared_ptr<json::Object> gridConfiguration) {
	// nullchecks
	// check json
	if (gridConfiguration == NULL) {
		glass3::util::Logger::log(
				"error",
				"CWeb::loadGridConfiguration: NULL json configuration.");
		return (false);
	}

	// grid definition variables and defaults
	std::string name = "UNKNOWN";
	int detect = CGlass::getNumStationsPerNode();
	int nucleate = CGlass::getNucleationDataCountThreshold();
	double thresh = CGlass::getNucleationStackThreshold();

	double resol = 0;
	double aziTaper = k_dAzimuthTaperDefault;
	double maxDepth = CGlass::k_dMaximumDepth;
	bool saveGrid = false;
	bool update = false;

	// get grid configuration from json
	// name
	if (((*gridConfiguration).HasKey("Name"))
			&& ((*gridConfiguration)["Name"].GetType()
					== json::ValueType::StringVal)) {
		name = (*gridConfiguration)["Name"].ToString();
	}

	// Nucleation Phases to be used for this generateLocalGrid
	if (((*gridConfiguration).HasKey("NucleationPhases"))
			&& ((*gridConfiguration)["NucleationPhases"].GetType()
					== json::ValueType::ObjectVal)) {
		json::Object nucleationPhases = (*gridConfiguration)["NucleationPhases"]
				.ToObject();

		if (loadTravelTimes(&nucleationPhases) == false) {
			glass3::util::Logger::log(
					"error",
					"CWeb::loadGridConfiguration: Error Loading NucleationPhases");
			return (false);
		}
	} else {
		if (loadTravelTimes((json::Object *) NULL) == false) {
			glass3::util::Logger::log(
					"error",
					"CWeb::loadGridConfiguration: Error Loading default "
					"NucleationPhases");
			return (false);
		}
	}

	// the number of stations that are assigned to each node
	if (((*gridConfiguration).HasKey("NumStationsPerNode"))
			&& ((*gridConfiguration)["NumStationsPerNode"].GetType()
					== json::ValueType::IntVal)) {
		detect = (*gridConfiguration)["NumStationsPerNode"].ToInt();
	}

	// number of picks that need to associate to start an event
	if (((*gridConfiguration).HasKey("NucleationDataCountThreshold"))
			&& ((*gridConfiguration)["NucleationDataCountThreshold"].GetType()
					== json::ValueType::IntVal)) {
		nucleate = (*gridConfiguration)["NucleationDataCountThreshold"].ToInt();
	}

	// viability threshold needed to exceed for a nucleation to be successful.
	if (((*gridConfiguration).HasKey("NucleationStackThreshold"))
			&& ((*gridConfiguration)["NucleationStackThreshold"].GetType()
					== json::ValueType::DoubleVal)) {
		thresh = (*gridConfiguration)["NucleationStackThreshold"].ToDouble();
	}

	// Node resolution for this grid
	if (((*gridConfiguration).HasKey("NodeResolution"))
			&& ((*gridConfiguration)["NodeResolution"].GetType()
					== json::ValueType::DoubleVal)) {
		resol = (*gridConfiguration)["NodeResolution"].ToDouble();
	} else {
		glass3::util::Logger::log(
				"error",
				"CWeb::loadGridConfiguration: Missing required Resolution Key.");
		return (false);
	}

	// sets the aziTaper value
	if ((*gridConfiguration).HasKey("AzimuthGapTaper")
			&& ((*gridConfiguration)["AzimuthGapTaper"].GetType()
					== json::ValueType::DoubleVal)) {
		aziTaper = (*gridConfiguration)["AzimuthGapTaper"].ToDouble();
	}

	// sets the maxDepth value
	if ((*gridConfiguration).HasKey("MaximumDepth")
			&& ((*gridConfiguration)["MaximumDepth"].GetType()
					== json::ValueType::DoubleVal)) {
		maxDepth = (*gridConfiguration)["MaximumDepth"].ToDouble();
	}

	// whether to create a file detailing the node configuration for
	// this grid
	if ((gridConfiguration->HasKey("SaveGrid"))
			&& ((*gridConfiguration)["SaveGrid"].GetType()
					== json::ValueType::BoolVal)) {
		saveGrid = (*gridConfiguration)["SaveGrid"].ToBool();
	}

	// set whether to update weblists
	if ((gridConfiguration->HasKey("UpdateGrid"))
			&& ((*gridConfiguration)["UpdateGrid"].GetType()
					== json::ValueType::BoolVal)) {
		update = (*gridConfiguration)["UpdateGrid"].ToBool();
	}

	// set whether to use a zonestats file
	if ((gridConfiguration->HasKey("ZoneStatsFile"))
			&& ((*gridConfiguration)["ZoneStatsFile"].GetType()
					== json::ValueType::StringVal)) {
		m_sZoneStatsFileName = (*gridConfiguration)["ZoneStatsFile"].ToString();
	} else {
		m_sZoneStatsFileName = "";
	}

	// sets the m_dDepthResolution value
	if ((*gridConfiguration).HasKey("DepthResolution")
			&& ((*gridConfiguration)["DepthResolution"].GetType()
					== json::ValueType::DoubleVal)) {
		m_dDepthResolution = (*gridConfiguration)["DepthResolution"].ToDouble();
	} else {
		m_dDepthResolution = k_dDepthResolutionUndefined;
	}

	// initialize
	initialize(name, thresh, detect, nucleate, resol, update, saveGrid,
				m_pNucleationTravelTime1, m_pNucleationTravelTime2, aziTaper,
				maxDepth);

	// generate site and network filter lists
	loadSiteFilters(gridConfiguration);

	// Generate eligible station list
	loadWebSiteList();

	// Load zone statistics
	if (!m_sZoneStatsFileName.empty()) {
		m_pZoneStats = std::make_shared<traveltime::CZoneStats>();
		if (!m_pZoneStats->setup(&m_sZoneStatsFileName)) {
			glass3::util::Logger::log(
					"error",
					"CWeb::loadGridConfiguration: ZoneStats filename specified, "
							"but unable to load valid data from"
							+ m_sZoneStatsFileName);
			return (false);
		}
	}  // end if m_sZoneStatsFileName

	return (true);
}

// ---------------------------------------------------------loadTravelTimes
bool CWeb::loadTravelTimes(json::Object *gridConfiguration) {
	std::lock_guard<std::recursive_mutex> ttGuard(m_travelTimeMutex);
	// check json
	if (gridConfiguration == NULL) {
		glass3::util::Logger::log(
				"info", "CWeb::loadTravelTimes: NULL json configuration "
				"Using default first phase and no second phase");

		// if no json object, default to P
		// clean out old phase if any
		m_pNucleationTravelTime1.reset();

		// use overall glass default if available
		if (CGlass::getDefaultNucleationTravelTime() != NULL) {
			m_pNucleationTravelTime1 = CGlass::getDefaultNucleationTravelTime();
		} else {
			// create new traveltime
			m_pNucleationTravelTime1 =
					std::make_shared<traveltime::CTravelTime>();

			// set up the nucleation traveltime
			m_pNucleationTravelTime1->setup("P");
		}

		// no second phase
		// clean out old phase if any
		m_pNucleationTravelTime2.reset();

		return (true);
	}

	// load the first travel time
	if ((gridConfiguration->HasKey("Phase1"))
			&& ((*gridConfiguration)["Phase1"].GetType()
					== json::ValueType::ObjectVal)) {
		// get the phase object
		json::Object phsObj = (*gridConfiguration)["Phase1"].ToObject();

		// clean out old phase if any
		m_pNucleationTravelTime1.reset();

		// create new traveltime
		m_pNucleationTravelTime1 = std::make_shared<traveltime::CTravelTime>();

		// get the phase name
		// default to P
		std::string phs = "P";
		if (phsObj.HasKey("PhaseName")) {
			phs = phsObj["PhaseName"].ToString();
		}

		// get the file if present
		std::string file = "";
		if (phsObj.HasKey("TravFile")) {
			file = phsObj["TravFile"].ToString();

			glass3::util::Logger::log(
					"info",
					"CWeb::loadTravelTimes: Using file location: " + file
							+ " for first phase: " + phs);
		} else {
			glass3::util::Logger::log(
					"info",
					"CWeb::loadTravelTimes: Using default file location for "
							"first phase: " + phs);
		}

		// set up the first phase travel time
		if (m_pNucleationTravelTime1->setup(phs, file) == false) {
			glass3::util::Logger::log(
					"error", "CWeb::loadTravelTimes: Failed to load file "
							"location " + file + " for first phase: " + phs);
			return (false);
		}

		std::string ttDebugPath = "./";
		if (phsObj.HasKey("TTDebugPath")) {
			ttDebugPath = phsObj["TTDebugPath"].ToString();
		}

		int numTTDebugDepths = 0;
		json::Array ttDebugDepths;
		if ((phsObj.HasKey("TTDebugDepth") &&
				(phsObj["TTDebugDepth"].GetType() == json::ValueType::ArrayVal))) {
			ttDebugDepths = phsObj["TTDebugDepth"].ToArray();
			numTTDebugDepths = ttDebugDepths.size();
		}

		if ((ttDebugPath != "") && (numTTDebugDepths > 0)) {
			for (int z = 0; z < numTTDebugDepths; z++) {
				std::string fileName = ttDebugPath + "/" + m_sName +
					+ "nuc_" + m_pNucleationTravelTime1->m_sPhase	+ "_"
					+ std::to_string(ttDebugDepths[z].ToDouble()) + ".csv";

				m_pNucleationTravelTime1->writeToFile(fileName,
					ttDebugDepths[z].ToDouble());
			}
		}
	} else {
		// if no first phase, use default from glass
		m_pNucleationTravelTime1.reset();

		// use overall glass default if available
		if (CGlass::getDefaultNucleationTravelTime() != NULL) {
			m_pNucleationTravelTime1 = CGlass::getDefaultNucleationTravelTime();
		} else {
			// create new traveltime
			m_pNucleationTravelTime1 =
					std::make_shared<traveltime::CTravelTime>();

			// set up the nucleation traveltime
			m_pNucleationTravelTime1->setup("P");
		}

		glass3::util::Logger::log(
				"info", "CWeb::loadTravelTimes: Using default first phase");
	}

	// load the second travel time
	if ((gridConfiguration->HasKey("Phase2"))
			&& ((*gridConfiguration)["Phase2"].GetType()
					== json::ValueType::ObjectVal)) {
		// get the phase object
		json::Object phsObj = (*gridConfiguration)["Phase2"].ToObject();

		// clean out old phase if any
		m_pNucleationTravelTime2.reset();

		// create new travel time
		m_pNucleationTravelTime2 = std::make_shared<traveltime::CTravelTime>();

		// get the phase name
		// default to S
		std::string phs = "S";
		if (phsObj.HasKey("PhaseName")) {
			phs = phsObj["PhaseName"].ToString();
		}

		// get the file if present
		std::string file = "";
		if (phsObj.HasKey("TravFile")) {
			file = phsObj["TravFile"].ToString();

			glass3::util::Logger::log(
					"info",
					"CWeb::loadTravelTimes: Using file location: " + file
							+ " for second phase: " + phs);
		} else {
			glass3::util::Logger::log(
					"info",
					"CWeb::loadTravelTimes: Using default file location for "
							"second phase: " + phs);
		}

		// set up the second phase travel time
		if (m_pNucleationTravelTime2->setup(phs, file) == false) {
			glass3::util::Logger::log(
					"error", "CWeb::loadTravelTimes: Failed to load file "
							"location " + file + " for second phase: " + phs);
			return (false);
		}

		std::string ttDebugPath = "./";
		if (phsObj.HasKey("TTDebugPath")) {
			ttDebugPath = phsObj["TTDebugPath"].ToString();
		}

		int numTTDebugDepths = 0;
		json::Array ttDebugDepths;
		if ((phsObj.HasKey("TTDebugDepth") &&
				(phsObj["TTDebugDepth"].GetType() == json::ValueType::ArrayVal))) {
			ttDebugDepths = phsObj["TTDebugDepth"].ToArray();
			numTTDebugDepths = ttDebugDepths.size();
		}

		if ((ttDebugPath != "") && (numTTDebugDepths > 0)) {
			for (int z = 0; z < numTTDebugDepths; z++) {
				std::string fileName = ttDebugPath + "/" + m_sName +
					+ "nuc_" + m_pNucleationTravelTime2->m_sPhase	+ "_"
					+ std::to_string(ttDebugDepths[z].ToDouble()) + ".csv";

				m_pNucleationTravelTime2->writeToFile(fileName,
					ttDebugDepths[z].ToDouble());
			}
		}
	} else {
		// no second phase
		// clean out old phase if any
		m_pNucleationTravelTime2.reset();

		glass3::util::Logger::log(
				"info",
				"CWeb::loadTravelTimes: Not using secondary nucleation phase");
	}

	return (true);
}

// ---------------------------------------------------------loadSiteFilters
bool CWeb::loadSiteFilters(std::shared_ptr<json::Object> gridConfiguration) {
	// nullchecks
	// check json
	if (gridConfiguration == NULL) {
		glass3::util::Logger::log("error",
									"genSiteFilters: NULL json configuration.");
		return (false);
	}

	std::lock_guard<std::recursive_mutex> webGuard(m_WebMutex);

	// Get the network names to be included in this web.
	if ((*gridConfiguration).HasKey("IncludeNetworks")
			&& ((*gridConfiguration)["IncludeNetworks"].GetType()
					== json::ValueType::ArrayVal)) {
		// clear any previous network filter list
		m_vNetworksFilter.clear();

		// get the network array
		json::Array arr = (*gridConfiguration)["IncludeNetworks"].ToArray();

		// for each network in the array
		int netFilterCount = 0;
		for (int i = 0; i < arr.size(); i++) {
			if (arr[i].GetType() == json::ValueType::StringVal) {
				// get the network
				std::string snet = arr[i].ToString();

				// add to the network filter list
				m_vNetworksFilter.push_back(snet);
				netFilterCount++;
			}
		}

		if (netFilterCount > 0) {
			glass3::util::Logger::log(
					"debug",
					"CWeb::genSiteFilters: " + std::to_string(netFilterCount)
							+ " network filters configured.");
		}
	}

	// get the SCNL names to be included in this web.
	if ((*gridConfiguration).HasKey("IncludeSites")
			&& ((*gridConfiguration)["IncludeSites"].GetType()
					== json::ValueType::ArrayVal)) {
		// clear any previous site filter list
		m_vSitesFilter.clear();

		// get the site array
		int staFilterCount = 0;
		json::Array arr = (*gridConfiguration)["IncludeSites"].ToArray();

		// for each site in the array
		for (int i = 0; i < arr.size(); i++) {
			if (arr[i].GetType() == json::ValueType::StringVal) {
				// get the scnl
				std::string scnl = arr[i].ToString();

				// add to the site filter list
				m_vSitesFilter.push_back(scnl);
				staFilterCount++;
			}
		}

		if (staFilterCount > 0) {
			glass3::util::Logger::log(
					"debug",
					"CWeb::genSiteFilters: " + std::to_string(staFilterCount)
							+ " SCNL filters configured.");
		}
	}

	// check to see if we're only to use teleseismic stations.
	if ((*gridConfiguration).HasKey("UseOnlyTeleseismicStations")
			&& ((*gridConfiguration)["UseOnlyTeleseismicStations"].GetType()
					== json::ValueType::BoolVal)) {
		m_bUseOnlyTeleseismicStations =
				(*gridConfiguration)["UseOnlyTeleseismicStations"].ToBool();

		glass3::util::Logger::log(
				"debug",
				"CWeb::genSiteFilters: bUseOnlyTeleseismicStations is "
						+ std::to_string(m_bUseOnlyTeleseismicStations) + ".");
	}

	// check to see if we're only to use teleseismic stations.
	if ((*gridConfiguration).HasKey("QualityFilter")
			&& ((*gridConfiguration)["QualityFilter"].GetType()
					== json::ValueType::DoubleVal)) {
		m_dQualityFilter =
				(*gridConfiguration)["QualityFilter"].ToDouble();

		glass3::util::Logger::log(
				"debug",
				"CWeb::genSiteFilters: m_dQualityFilter is "
						+ std::to_string(m_dQualityFilter) + ".");
	} else {
		m_dQualityFilter = -1.0;
	}

	return (true);
}

// ---------------------------------------------------------isSiteAllowed
bool CWeb::isSiteAllowed(std::shared_ptr<CSite> site, bool checkEnabled) {
	if (site == NULL) {
		return (false);
	}

	if (checkEnabled == true) {
		// if the site has been disabled, return false
		if (!site->getUse()) {
			return(false);
		}
		if (!site->getEnable()) {
			return(false);
		}
		// check quality
		if (m_dQualityFilter > 0) {
			if (site->getQuality() < m_dQualityFilter) {
				return(false);
			}
		}
	}

	// If we have do not have a site and/or network filter, just add it
	if ((m_vNetworksFilter.size() == 0) && (m_vSitesFilter.size() == 0)
			&& (m_bUseOnlyTeleseismicStations == false)) {
		return (true);
	}

	// default to false
	bool returnFlag = false;

	// if we have a network filter, make sure network is allowed before adding
	if ((m_vNetworksFilter.size() > 0)
			&& (find(m_vNetworksFilter.begin(), m_vNetworksFilter.end(),
						site->getNetwork()) != m_vNetworksFilter.end())) {
		returnFlag = true;
	}

	// if we have a site filter, make sure site is allowed before adding
	if ((m_vSitesFilter.size() > 0)
			&& (find(m_vSitesFilter.begin(), m_vSitesFilter.end(),
						site->getSCNL()) != m_vSitesFilter.end())) {
		returnFlag = true;
	}

	// if we're only using teleseismic stations, make sure site is allowed
	// before adding
	if (m_bUseOnlyTeleseismicStations == true) {
		// is this site used for teleseismic
		if (site->getUseForTeleseismic() == true) {
			returnFlag = true;
		} else {
			returnFlag = false;
		}
	}

	return (returnFlag);
}

// ---------------------------------------------------------loadWebSiteList
bool CWeb::loadWebSiteList() {
	// nullchecks
	// check pSiteList
	if (m_pSiteList == NULL) {
		glass3::util::Logger::log(
				"error", "CWeb::loadWebSiteList: NULL pSiteList pointer.");
		return (false);
	}

	// check to see if the site list has changed since we last asked
	// this is to reduce the number of redundant site list updates
	int tUpdated = m_pSiteList->getLastUpdated();
	if (tUpdated > m_tLastUpdated) {
		m_tLastUpdated = tUpdated;
	} else {
		// our site list is already the latest
		return(true);
	}

	// git all the sites in pSiteList
	std::vector<std::shared_ptr<CSite>> siteList =
			m_pSiteList->getListOfSites();

	// get the total number sites in glass's site list
	int nsite = siteList.size();

	// don't bother continuing if we have no sites
	if (nsite <= 0) {
		glass3::util::Logger::log(
				"warning", "CWeb::loadWebSiteList: No sites in site list.");
		return (false);
	}

	// log
	// char sLog[glass3::util::Logger::k_nMaxLogEntrySize];
	// snprintf(sLog, sizeof(sLog),
	// "CWeb::loadWebSiteList: %d sites available for web %s", nsite,
	// m_sName.c_str());
	// glass3::util::Logger::log("debug", sLog);

	// lock the site list while we are rebuilding it
	while ((m_vSiteMutex.try_lock() == false) &&
					(getTerminate() == false)) {
		// update thread status
		setThreadHealth(true);
	}

	// clear web site list
	vWebSites.clear();

	// update thread status
	setThreadHealth(true);

	// for each site
	for (std::vector<std::shared_ptr<CSite>>::iterator it = siteList.begin();
			it != siteList.end(); ++it) {
		// get site from the overall site list
		std::shared_ptr<CSite> site = *it;

		// nullcheck
		if (site == NULL) {
			continue;
		}

		// if this site meets our criteria, add it to the sorted list
		if (isSiteAllowed(site)) {
			vWebSites.push_back(
					std::pair<double, std::shared_ptr<CSite>>(0.0, site));
		}
	}

	// log
	// snprintf(sLog, sizeof(sLog),
	// "CWeb::loadWebSiteList: selected %d allowed sites for web %s",
	// vWebSites.size(), m_sName.c_str());
	// glass3::util::Logger::log("debug", sLog);

	m_vSiteMutex.unlock();

	return (true);
}

// ---------------------------------------------------------sortSiteListForNode
std::vector<std::pair<double, std::shared_ptr<CSite>>>
	CWeb::sortSiteListForNode(double lat, double lon, double depth,
		std::vector<std::pair<double, std::shared_ptr<CSite>>> sites) {
	// set to provided geographic location
	glass3::util::Geo geo;
	geo.setGeographic(lat, lon, glass3::util::Geo::k_EarthRadiusKm - depth);

	// compute the distance to each site
	for (std::vector<std::pair<double, std::shared_ptr<CSite>>>::iterator it
		= sites.begin();
		it != sites.end(); ++it) {
		// get the site distance pair
		auto sitePair = *it;

		// get the site from the second part of the pair
		std::shared_ptr<CSite> site = sitePair.second;

		// nullcheck
		if (site == NULL) {
			continue;
		}

		// compute the distance
		sitePair.first = site->getDelta(&geo);

		// set the distance in the vector
		*it = sitePair;
	}

	// update thread status
	setThreadHealth(true);

	// sort site distance pair list
	std::sort(sites.begin(), sites.end(), sortSite);

	return(sites);
}

// ---------------------------------------------------------generateNode
std::shared_ptr<CNode> CWeb::generateNode(double lat, double lon, double z,
		double resol, std::vector<std::pair<double, std::shared_ptr<CSite>>> sites) {
	// nullcheck
	if ((m_pNucleationTravelTime1 == NULL)
			&& (m_pNucleationTravelTime2 == NULL)) {
		glass3::util::Logger::log("error",
									"CWeb::genNode: No valid trav pointers.");
		return (NULL);
	}

	double maxZ = m_dMaxDepth;
	if (m_dDepthResolution != -1) {
		maxZ = z + m_dDepthResolution;
	}

	// create node
	std::shared_ptr<CNode> node(new CNode(m_sName, lat, lon, z, resol, maxZ));

	// set parent web
	node->setWeb(this);

	// return empty node if we don't
	// have any sites
	if (sites.size() == 0) {
		return (node);
	}

	// generate the sites for the node
	node = generateNodeSites(node, sites);

	// return the populated node
	return (node);
}

// ---------------------------------------------------------addNode
bool CWeb::addNode(std::shared_ptr<CNode> node) {
	// nullcheck
	if (node == NULL) {
		return (false);
	}
	std::lock_guard<std::mutex> vNodeGuard(m_vNodeMutex);

	m_vNode.push_back(node);

	return (true);
}

// ---------------------------------------------------------generateNodeSites
std::shared_ptr<CNode> CWeb::generateNodeSites(std::shared_ptr<CNode> node,
	std::vector<std::pair<double, std::shared_ptr<CSite>>> sites) {
	// nullchecks
	// check node
	if (node == NULL) {
		glass3::util::Logger::log("error",
									"CWeb::genNodeSites: NULL node pointer.");
		return (NULL);
	}
	// check trav
	if ((m_pNucleationTravelTime1 == NULL)
			&& (m_pNucleationTravelTime2 == NULL)) {
		glass3::util::Logger::log(
				"error", "CWeb::genNodeSites: No valid trav pointers.");
		return (NULL);
	}
	// check sites
	if (sites.size() == 0) {
		glass3::util::Logger::log("error", "CWeb::genNodeSites: No sites.");
		return (node);
	}
	// check nDetect
	if (m_iNumStationsPerNode == 0) {
		glass3::util::Logger::log("error", "CWeb::genNodeSites: num "
			"stations per node is 0.");
		return (node);
	}

	// update thread status
	setThreadHealth(true);

	int sitesAllowed = m_iNumStationsPerNode;
	int sitesAvailable = sites.size();
	if (sitesAvailable < sitesAllowed) {
		glass3::util::Logger::log("warning",
									"CWeb::genNodeSites: the number of sites allowed ("
									+ std::to_string(sitesAllowed) +
									") is greater than the number of sites available (" +
									std::to_string(sitesAvailable) + ") for web "
									+ m_sName + " using the number of sites available.");
		sitesAllowed = sitesAvailable;
	}

	// we're now modifying a node
	node->setEnabled(false);

	// clear node of any existing sites (on adds/removes/updates)
	node->clearSiteLinks();

	// for all the allowed sites (note we won't usually end up going through all
	// of them)
	int siteCount = 0;
	for (std::vector<std::pair<double, std::shared_ptr<CSite>>>::iterator it
			= sites.begin();
			it != sites.end(); ++it) {
		// break out of loop if we've exceeded
		// the number of allowed sites
		if (siteCount > sitesAllowed) {
			break;
		}

		// get the site distance pair
		auto sitePair = *it;

		// get the site from the second part of the pair
		std::shared_ptr<CSite> site = sitePair.second;

		// nullcheck site
		if (site == NULL) {
			// move on to the next site
			continue;
		}

		// add current site to this node
		if (addSiteToNode(site, node) == false) {
			// move on to the next site
			continue;
		} else {
			// we've added a site
			siteCount++;
		}
	}

	// sort the site links in ascending distance
	node->sortSiteLinks();

	// done with node
	node->setEnabled(true);

	// return updated node
	return (node);
}

// ---------------------------------------------------------addSiteToNode
bool CWeb::addSiteToNode(std::shared_ptr<CSite> newSite,
	std::shared_ptr<CNode> node) {
	if (newSite == NULL) {
		return(false);
	}

	if (node == NULL) {
		return(false);
	}

	// update thread status
	setThreadHealth(true);

	// compute distance between site and node
	glass3::util::Geo geoNode = node->getGeo();
	double distance = glass3::util::GlassMath::k_RadiansToDegrees
			* newSite->getDelta(&geoNode);

	// check distance
	if (distance < 0) {
		glass3::util::Logger::log(
				"error",
				"CWeb::addSiteToNode: Invalid distance for site "
						+ newSite->getSCNL() + " for node "
						+ std::to_string(node->getLatitude())
						+ "," + std::to_string(node->getLongitude())
						+ "," + std::to_string(node->getDepth())
						+ " in web " + m_sName + ".");

		return(false);
	}

	// lock the travel times so it is constant while we build the node
	while ((m_travelTimeMutex.try_lock() == false) &&
				 (getTerminate() == false)) {
		// update thread status
		setThreadHealth(true);

		// wait a little while
		std::this_thread::sleep_for(
				std::chrono::milliseconds(getSleepTime()));
	}

	// setup traveltimes for this node
	if (m_pNucleationTravelTime1 != NULL) {
		m_pNucleationTravelTime1->setTTOrigin(node->getLatitude(),
											node->getLongitude(),
											node->getDepth());
	}
	if (m_pNucleationTravelTime2 != NULL) {
		m_pNucleationTravelTime2->setTTOrigin(node->getLatitude(),
											node->getLongitude(),
											node->getDepth());
	}

	// compute traveltime(s) between site and node using distance
	// we assume that the caller initialized the travel times for us
	double travelTime1 = traveltime::CTravelTime::k_dTravelTimeInvalid;
	std::string phase1 = traveltime::CTravelTime::k_dPhaseInvalid;
	if (m_pNucleationTravelTime1 != NULL) {
		travelTime1 = m_pNucleationTravelTime1->T(distance);
		phase1 = m_pNucleationTravelTime1->m_sPhase;
	}

	double travelTime2 = traveltime::CTravelTime::k_dTravelTimeInvalid;
	std::string phase2 = traveltime::CTravelTime::k_dPhaseInvalid;
	if (m_pNucleationTravelTime2 != NULL) {
		travelTime2 = m_pNucleationTravelTime2->T(distance);
		phase2 = m_pNucleationTravelTime2->m_sPhase;
	}

	m_travelTimeMutex.unlock();

	// check to make sure we have at least one valid travel time
	if ((travelTime1 < 0) && (travelTime2 < 0)) {
		glass3::util::Logger::log(
				"warning",
				"CWeb::addSiteToNode: No valid travel times for site "
						+ newSite->getSCNL() + " for " + phase1 + " or " + phase2
						+ " at distance: " + std::to_string(distance)
						+ " for node " + std::to_string(node->getLatitude())
						+ "," + std::to_string(node->getLongitude())
						+ "," + std::to_string(node->getDepth())
						+ " in web " + m_sName + ".");

		return(false);
	}

	// Link node to new site using traveltime(s)
	if (node->linkSite(newSite, node, distance, travelTime1, phase1,
						travelTime2, phase2) == false) {
		glass3::util::Logger::log(
				"error",
				"CWeb::addSiteToNode: Failed to link site "
						+ newSite->getSCNL() + " to node "
						+ std::to_string(node->getLatitude())
						+ "," + std::to_string(node->getLongitude())
						+ " for web " + m_sName + ".");
		return(false);
	}

	// update thread status
	setThreadHealth(true);

	return(true);
}

// ---------------------------------------------------------updateSite
void CWeb::updateSite(std::shared_ptr<CSite> site) {
	//  nullcheck
	if (site == NULL) {
		return;
	}

	// don't bother if we're not allowed to update
	if (m_bUpdate == false) {
		return;
	}

	// don't bother if this site isn't allowed, note that we ignore the
	// use / enabled /etc flags here since we're updating a site (which might be
	// unused or disabled, but still in this web currently)
	if (isSiteAllowed(site, false) == false) {
		glass3::util::Logger::log(
				"debug",
				"CWeb::updateSite: Station " + site->getSCNL()
						+ " not allowed in web " + m_sName + ".");
		return;
	}

	// log what we're doing
	if (site->getUse() == true) {
		glass3::util::Logger::log(
			"debug",
			"CWeb::updateSite: Trying to add/update site " + site->getSCNL()
					+ " in web: " + m_sName + ".");
	} else {
		glass3::util::Logger::log(
			"debug",
			"CWeb::updateSite: Trying to remove site " + site->getSCNL()
					+ " from web: " + m_sName + ".");
	}

	int nodeModCount = 0;
	int nodeCount = 0;
	int totalNodes = size();

	std::chrono::high_resolution_clock::time_point tStartTime =
			std::chrono::high_resolution_clock::now();

	// lock the site list while getting site list
	while ((m_vSiteMutex.try_lock() == false) &&
					(getTerminate() == false)) {
		// update thread status
		setThreadHealth(true);

		// wait a little while
		std::this_thread::sleep_for(
				std::chrono::milliseconds(getSleepTime()));
	}

	// update the web site list from the master list in SiteList
	loadWebSiteList();

	// grab a copy of the site vector for our use
	std::vector<std::pair<double, std::shared_ptr<CSite>>>
		sites(vWebSites);

	m_vSiteMutex.unlock();

	// for each node in web
	for (auto &node : m_vNode) {
		nodeCount++;

		// update thread status
		setThreadHealth(true);

		// some optimization to avoid unneeded work in some cases
		std::string whatDo = "";
		if (site->getUse() == true) {
			// use/enable true means add/update
			// do we have this site
			std::shared_ptr<CSite> oldSite = node->getSite(site->getSCNL());
			if (oldSite == NULL) {
				// we don't have it, so it's an add
				whatDo = "Add";

				// if we're an add check to make sure the new station is valid
				// for the node i.e. not too far out,
				glass3::util::Geo geo;
				geo.setGeographic(node->getLatitude(), node->getLongitude(),
									glass3::util::Geo::k_EarthRadiusKm - node->getDepth());

				// compute delta distance between site and node
				double newDistance = glass3::util::GlassMath::k_RadiansToDegrees
						* site->getGeo().delta(&geo);

				// get last site in node list
				// NOTE: this assumes that the node site list is sorted
				// on distance
				std::shared_ptr<CSite> furthestSite = node->getLastSite();

				// compute distance to farthest site
				double maxDistance = glass3::util::GlassMath::k_RadiansToDegrees
						* geo.delta(&furthestSite->getGeo());

				// Ignore if new site is farther than last linked site, and
				// if we're at max sites (no more space)
				if ((node->getSiteLinksCount() >= m_iNumStationsPerNode)
						&& (newDistance > maxDistance)) {
					continue;
				}
			} else {
				// we already have it, so it's an update
				whatDo = "Update";

				// now check again if the site is allowed,
				// INCLUDING use / enabled / etc this time
				if (isSiteAllowed(site) == true) {
					// if it still belongs in the web, then there's nothing to do
					// we may have already done this update with a previous modification
					// NOTE we trust that the site hasn't MOVED...
					continue;
				}
				// can't check further if sites are the same, because sitelist has already
				// updated the site, so it'll always appear the same at this point,
				// so all we can do is rebuild the node in case the updated
				// changes are important to the node.
			}
		} else {
			// use/enable false means remove
			whatDo = "Remove";

			// don't bother if this node doesn't have this site, there's
			// nothing to remove here
			if (node->getSite(site->getSCNL()) == NULL) {
				continue;
			}
		}

		// resort the site list for the current node,
		sites = sortSiteListForNode(node->getLatitude(), node->getLongitude(),
			node->getDepth(), sites);

		// its easier to just regenerate the node links than rat through
		// and figure out where to add / remove / update the new station.
		// so rebuild the node
		generateNodeSites(node, sites);

		// we've modified a node
		nodeModCount++;

		// modding by 1000 ensures we don't get that many log entries
		if (nodeCount % 1000 == 0) {
			glass3::util::Logger::log(
					"debug",
					"CWeb::updateSite: Station " + site->getSCNL() + " processed "
							+ "(" + whatDo + ") "
							+ std::to_string(nodeCount) + " out of "
							+ std::to_string(totalNodes) + " nodes in web: "
							+ m_sName + ". Modified "
							+ std::to_string(nodeModCount) + " nodes so far.");
		}
	}

	std::chrono::high_resolution_clock::time_point tEndTime =
			std::chrono::high_resolution_clock::now();

	double updateTime =
			std::chrono::duration_cast<std::chrono::duration<double>>(
					tEndTime - tStartTime).count();

	// log info if we've updated any node with this site
	if (nodeModCount > 0) {
		char sLog[glass3::util::Logger::k_nMaxLogEntrySize];
		snprintf(sLog, sizeof(sLog), "CWeb::updateSite: Station %s modified"
					" %d node(s) in web: %s in %f seconds",
					site->getSCNL().c_str(), nodeModCount, m_sName.c_str(), updateTime);
		glass3::util::Logger::log("info", sLog);
	} else {
		glass3::util::Logger::log(
				"debug",
				"CWeb::updateSite: Station " + site->getSCNL()
						+ " did not modify in any "
						"nodes in web: " + m_sName + " in " + std::to_string(updateTime)
						+ " seconds.");
	}
}

// ---------------------------------------------------------addJob
void CWeb::addJob(std::function<void()> newjob) {
	if (getNumThreads() == 0) {
		// no background thread, just run the job
		try {
			newjob();
		} catch (const std::exception &e) {
			glass3::util::Logger::log(
					"error",
					"CWeb::jobLoop: Exception during job(): "
							+ std::string(e.what()));
		}

		return;
	}

	// add the job to the queue
	std::lock_guard<std::mutex> guard(m_QueueMutex);
	m_JobQueue.push(newjob);
}

// ---------------------------------------------------------work
glass3::util::WorkState CWeb::work() {
	// lock for queue access
	m_QueueMutex.lock();

	// are there any jobs
	if (m_JobQueue.empty() == true) {
		// unlock and skip until next time
		m_QueueMutex.unlock();

		return (glass3::util::WorkState::Idle);
	}

	// get the next job
	std::function<void()> newjob = m_JobQueue.front();
	m_JobQueue.pop();

	// done with queue
	m_QueueMutex.unlock();

	// run the job
	try {
		newjob();
	} catch (const std::exception &e) {
		glass3::util::Logger::log(
				"error",
				"CWeb::workLoop: Exception during job(): "
						+ std::string(e.what()));
		return (glass3::util::WorkState::Error);
	}

	// give up some time at the end of the loop
	return (glass3::util::WorkState::OK);
}

// ---------------------------------------------------------hasSite
bool CWeb::hasSite(std::shared_ptr<CSite> site) {
	//  nullcheck
	if (site == NULL) {
		return (false);
	}

	std::lock_guard<std::mutex> vNodeGuard(m_vNodeMutex);

	// for each node in web
	for (auto &node : m_vNode) {
		// check to see if we have this site
		if (node->getSite(site->getSCNL()) != NULL) {
			// found it
			return (true);
		}
	}

	// site not found
	return (false);
}

// ---------------------------------------------------------getAzimuthTaper
double CWeb::getAzimuthTaper() const {
	return (m_dAzimuthTaper);
}

// ---------------------------------------------------------getMaxDepth
double CWeb::getMaxDepth() const {
	return (m_dMaxDepth);
}

// ---------------------------------------------------------getSiteList
const CSiteList* CWeb::getSiteList() const {
	std::lock_guard<std::recursive_mutex> webGuard(m_WebMutex);
	return (m_pSiteList);
}

// ---------------------------------------------------------setSiteList
void CWeb::setSiteList(CSiteList* siteList) {
	std::lock_guard<std::recursive_mutex> webGuard(m_WebMutex);
	m_pSiteList = siteList;
}

// ---------------------------------------------------------getUpdate
bool CWeb::getUpdate() const {
	return (m_bUpdate);
}

// ---------------------------------------------------------getSaveGrid
bool CWeb::getSaveGrid() const {
	return (m_bSaveGrid);
}

// ---------------------------------------------------------getResolution
double CWeb::getNodeResolution() const {
	return (m_dNodeResolution);
}

// --------------------------------------------------getNucleationStackThreshold
double CWeb::getNucleationStackThreshold() const {
	return (m_dNucleationStackThreshold);
}

// -------------------------------------------------------getNumStationsPerNode
int CWeb::getNumStationsPerNode() const {
	return (m_iNumStationsPerNode);
}

// ---------------------------------------------------getNucleationDataThreshold
int CWeb::getNucleationDataCountThreshold() const {
	return (m_iNucleationDataCountThreshold);
}

// ---------------------------------------------------------getName
const std::string& CWeb::getName() const {
	return (m_sName);
}

// -----------------------------------------------------getNucleationTravelTime1
const std::shared_ptr<traveltime::CTravelTime>& CWeb::getNucleationTravelTime1() const {  // NOLINT
	std::lock_guard<std::recursive_mutex> webGuard(m_WebMutex);
	return (m_pNucleationTravelTime1);
}

// -----------------------------------------------------getNucleationTravelTime2
const std::shared_ptr<traveltime::CTravelTime>& CWeb::getNucleationTravelTime2() const {  // NOLINT
	std::lock_guard<std::recursive_mutex> webGuard(m_WebMutex);
	return (m_pNucleationTravelTime2);
}

// -----------------------------------------------------getNetworksFilterSize
int CWeb::getNetworksFilterSize() const {
	std::lock_guard<std::recursive_mutex> webGuard(m_WebMutex);
	return (m_vNetworksFilter.size());
}

// ------------------------------------------------getUseOnlyTeleseismicStations
bool CWeb::getUseOnlyTeleseismicStations() const {
	std::lock_guard<std::recursive_mutex> webGuard(m_WebMutex);
	return (m_bUseOnlyTeleseismicStations);
}

// -----------------------------------------------------getSitesFilterSize
int CWeb::getSitesFilterSize() const {
	std::lock_guard<std::recursive_mutex> webGuard(m_WebMutex);
	return (m_vSitesFilter.size());
}

// -----------------------------------------------------size
int CWeb::size() const {
	std::lock_guard<std::mutex> vNodeGuard(m_vNodeMutex);
	return (m_vNode.size());
}
}  // namespace glasscore

