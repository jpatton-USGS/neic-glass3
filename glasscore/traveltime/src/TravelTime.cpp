#include "TravelTime.h"
#include <geo.h>
#include <logger.h>
#include <string>
#include <cmath>
#include <cstdio>
#include <cstring>
#include "TimeWarp.h"

namespace traveltime {

// constants
constexpr double CTravelTime::k_dTravelTimeInvalid;
const std::string CTravelTime::k_dPhaseInvalid = ""; // NOLINT

// structures
// Note: This structure matches the struct used to generate the travel time
// files and is needed to correctly read them
typedef struct _TTNucEntry
{
  float dTPhase; /* Time for this entry in seconds */
  float dDPhase; /* Distance for this entry in seconds */
  float dtdz; /* Vertical derivitive */
  float dtdx; /* Ray parameter (horizontal derivitive) */
  float dTOA; /* Take off angle */
  float dAssocStrength; /* binding strength to an origin */
  float dLocWeight; /* strength in affecting origin location */
  float dMagThresh; /* minimum Mx magnitude threshhold (not seen on smaller events) */ // NOLINT
  float dResidWidth; /* width of residual std-dev */
  char  * szPhase; /* phase name for this entry */
} TTNucEntry;

// the expected file format version
#define TT_FILE_FORMAT_VERSION 1003
#define TT_PHASE_NAME_SIZE 16

// ---------------------------------------------------------CTravelTime
CTravelTime::CTravelTime() {
	m_pDepthTravelTimeArray = NULL;
	m_pDepthDistanceArray = NULL;

	clear();
}

// ---------------------------------------------------------CTravelTime
CTravelTime::CTravelTime(const CTravelTime &travelTime) {
	m_pDepthTravelTimeArray = NULL;
	m_pDepthDistanceArray = NULL;

	clear();

	m_dDepth = travelTime.m_dDepth;
	m_dDelta = travelTime.m_dDepth;

	m_iNumDistances = travelTime.m_iNumDistances;
	m_dMinimumDistance = travelTime.m_dMinimumDistance;
	m_dMaximumDistance = travelTime.m_dMaximumDistance;

	m_iNumDepths = travelTime.m_iNumDepths;
	m_dMinimumDepth = travelTime.m_dMinimumDepth;
	m_dMaximumDepth = travelTime.m_dMaximumDepth;

	m_iNumTimes = travelTime.m_iNumTimes;
	m_dMinimumTime = travelTime.m_dMinimumTime;
	m_dMaximumTime = travelTime.m_dMaximumTime;

	m_dDistanceDelta = travelTime.m_dDistanceDelta;
	m_dDepthDelta = travelTime.m_dDepthDelta;
	m_dTimeDelta = travelTime.m_dTimeDelta;

	m_sPhase = travelTime.m_sPhase;

	m_pDepthDistanceArray = new double[m_iNumDistances * m_iNumDepths];
	m_pDepthTravelTimeArray = new double[m_iNumTimes * m_iNumDepths];

	for (int i = 0; i < (m_iNumDistances * m_iNumDepths); i++) {
		m_pDepthDistanceArray[i] = travelTime.m_pDepthDistanceArray[i];
	}

	for (int i = 0; i < (m_iNumTimes * m_iNumDepths); i++) {
		m_pDepthTravelTimeArray[i] = travelTime.m_pDepthTravelTimeArray[i];
	}
}

// ---------------------------------------------------------~CTravelTime
CTravelTime::~CTravelTime() {
	clear();
}

// ---------------------------------------------------------clear
void CTravelTime::clear() {
	m_dDepth = 0;
	m_dDelta = 0;

	m_iNumDistances = 0;
	m_dMinimumDistance = 0;
	m_dMaximumDistance = 0;

	m_iNumDepths = 0;
	m_dMinimumDepth = 0;
	m_dMaximumDepth = 0;

	m_iNumTimes = 0;
	m_dMinimumTime = 0;
	m_dMaximumTime = 0;

	m_dDistanceDelta = 0;
	m_dDepthDelta = 0;
	m_dTimeDelta = 0;

	m_sPhase = CTravelTime::k_dPhaseInvalid;

	if (m_pDepthTravelTimeArray) {
		delete (m_pDepthTravelTimeArray);
	}
	m_pDepthTravelTimeArray = NULL;

	if (m_pDepthDistanceArray) {
		delete (m_pDepthDistanceArray);
	}
	m_pDepthDistanceArray = NULL;
}


// ---------------------------------------------------------Setup
bool CTravelTime::setup(std::string phase, std::string file) {
	// nullcheck
	if (phase == CTravelTime::k_dPhaseInvalid) {
		glass3::util::Logger::log("error",
									"CTravelTime::Setup: empty phase provided");
		return (false);
	} else {
		m_sPhase = phase;
	}

	// generate file name if not specified
	if (file == "") {
		file = phase + ".ttt";
	}

	glass3::util::Logger::log(
			"debug", "CTravelTime::Setup: phase:" + phase + " file: " + file);

	// open file
	FILE *inFile = fopen(file.c_str(), "rb");
	if (!inFile) {
		glass3::util::Logger::log(
				"error", "CTravelTime::Setup: Cannot open file:" + file);
		return (false);
	}

	// first read the format version number
	int formatVersion = 0;
	fread(&formatVersion, sizeof(formatVersion), 1, inFile);

	if (formatVersion != TT_FILE_FORMAT_VERSION) {
		glass3::util::Logger::log(
				"error", "CTravelTime::Setup: Incorrect file format version: "
					+ std::to_string(formatVersion) + "for file: " + file);

		fclose(inFile);

		return (false);
	}

	// read the phase name
	char phaseName[TT_PHASE_NAME_SIZE];
	fread(phaseName, sizeof(phaseName), 1, inFile);
	m_sPhase = std::string(phaseName);

	// read the depth/distance attributes
	// read the number of distance points
	int numDist = 0;
	fread(&numDist, sizeof(numDist), 1, inFile);
	m_iNumDistances = numDist;

	// read the number of depth points
	int numDepth = 0;
	fread(&numDepth, sizeof(numDepth), 1, inFile);
	m_iNumDepths = numDepth;

	// read the minimum distance
	float minDist = 0;
	fread(&minDist, sizeof(minDist), 1, inFile);
	m_dMinimumDistance = static_cast<double>(minDist);

	// read the maximum distance
	float maxDist = 0;
	fread(&maxDist, sizeof(maxDist), 1, inFile);
	m_dMaximumDistance = static_cast<double>(maxDist);

	// read the minimum depth
	float minDepth = 0;
	fread(&minDepth, sizeof(minDepth), 1, inFile);
	m_dMinimumDepth = static_cast<double>(minDepth);

	// read the maximum depth
	float maxDepth = 0;
	fread(&maxDepth, sizeof(maxDepth), 1, inFile);
	m_dMaximumDepth = static_cast<double>(maxDepth);

	// read the distance delta
	float distDelta = 0;
	fread(&distDelta, sizeof(distDelta), 1, inFile);
	m_dDistanceDelta = static_cast<double>(distDelta);

	// read the depth delta
	float depthDelta = 0;
	fread(&depthDelta, sizeof(depthDelta), 1, inFile);
	m_dDepthDelta = static_cast<double>(depthDelta);

	// read the depth distance array
	// allocate the Depth/Distance table using the structure
  TTNucEntry * depthDistanceTable = reinterpret_cast<TTNucEntry *>(malloc(
		m_iNumDistances	* m_iNumDepths * sizeof(TTNucEntry)));

	// check the allocation
  if (!depthDistanceTable) {
    glass3::util::Logger::log(
				"error", "CTravelTime::Setup: ERROR! could not alloc %d bytes for "
				"Depth/Distance table with %d rows and %d cols from file(%s).",
                 (m_iNumDistances	* m_iNumDepths * sizeof(TTNucEntry)),
								 m_iNumDepths, m_iNumDistances, file);
		return(false);
  }

  // read the Depth/Distance table using the structure
	int numRead = fread(depthDistanceTable, sizeof(TTNucEntry),
		m_iNumDistances	* m_iNumDepths, inFile);

	// check how much we read
  if (numRead != (m_iNumDistances * m_iNumDepths)) {
    glass3::util::Logger::log(
				"error", "CTravelTime::Setup: ERROR! could not read entire "
				"Depth/Distance table from file(%s).  Expecting (%d) entries, read(%d)",
                file, (m_iNumDistances * m_iNumDepths), numRead);
		return(false);
  }

	// convert (copy out the relevant info from) the Depth/Distance table
	m_pDepthDistanceArray = new double[m_iNumDistances * m_iNumDepths];
	for (int index = 0; index < (m_iNumDistances * m_iNumDepths); index++) {
		m_pDepthDistanceArray[index] =
			static_cast<double>(depthDistanceTable[index].dDPhase);
	}

	// read the depth/time attributes
	// (we read the depth attributes above)
	// read the number of time points
	int numTime = 0;
	fread(&numDist, sizeof(numTime), 1, inFile);
	m_iNumTimes = numTime;

	// read the minimum time
	float minTime = 0;
	fread(&minTime, sizeof(minTime), 1, inFile);
	m_dMinimumTime = static_cast<double>(minTime);

	// read the maximum time
	float maxTime = 0;
	fread(&maxTime, sizeof(maxTime), 1, inFile);
	m_dMaximumTime = static_cast<double>(maxTime);

	// calculate the time delta
	m_dTimeDelta = (m_dMaximumTime - m_dMinimumTime) / (m_iNumTimes - 1);

	// read the depth time array
	// allocate the Depth/Time table using the structure
  TTNucEntry * depthTimeTable = reinterpret_cast<TTNucEntry *>(malloc(
		m_iNumTimes	* m_iNumDepths * sizeof(TTNucEntry)));

	// check the allocation
  if (!depthTimeTable) {
    glass3::util::Logger::log(
				"error", "CTravelTime::Setup: ERROR! could not alloc %d bytes for "
				"Depth/Time table with %d rows and %d cols from file(%s).",
                 (m_iNumTimes	* m_iNumDepths * sizeof(TTNucEntry)),
								 m_iNumDepths, m_iNumTimes, file);
		return(false);
  }

  // read the Depth/Time table using the structure
	int numRead = fread(depthTimeTable, sizeof(TTNucEntry),
		m_iNumTimes	* m_iNumDepths, inFile);

	// check how much we read
  if (numRead != (m_iNumTimes * m_iNumDepths)) {
    glass3::util::Logger::log(
				"error", "CTravelTime::Setup: ERROR! could not read entire "
				"Depth/Time table from file(%s).  Expecting (%d) entries, read(%d)",
                file, (m_iNumTimes * m_iNumDepths), numRead);
		return(false);
  }

	// convert (copy out the relevant info from) the Depth/Time table
	m_pDepthTravelTimeArray = new double[m_iNumTimes * m_iNumDepths];
	for (int index = 0; index < (m_iNumTimes * m_iNumDepths); index++) {
		m_pDepthTravelTimeArray[index] =
			static_cast<double>(depthTimeTable[index].dTPhase);
	}

	// finally read the format version number again to guard against
	// file corruption
	int endFormatVersion = 0;
	fread(&endFormatVersion, sizeof(endFormatVersion), 1, inFile);

	if (endFormatVersion != TT_FILE_FORMAT_VERSION) {
		glass3::util::Logger::log(
				"error", "CTravelTime::Setup: Incorrect end file format version: "
					+ std::to_string(formatVersion) + "for file: " + file
					+ ", possible corruption.");

		fclose(inFile);

		return (false);
	}

	// done with the file.
	fclose(inFile);

	return (true);
}

// ---------------------------------------------------------setOrigin
void CTravelTime::setTTOrigin(double lat, double lon, double depth) {
	m_geoTTOrigin.setGeographic(lat, lon,
								glass3::util::Geo::k_EarthRadiusKm - depth);
	m_dDepth = depth;
}

// ---------------------------------------------------------setOrigin
void CTravelTime::setTTOrigin(const glass3::util::Geo &geoOrigin) {
	m_geoTTOrigin = geoOrigin;
	// ditch dDepth in favor or
	m_dDepth = glass3::util::Geo::k_EarthRadiusKm
			- m_geoTTOrigin.m_dGeocentricRadius;
}

// ---------------------------------------------------------T
double CTravelTime::T(glass3::util::Geo *geo) {
	// Calculate travel time given CGeo
	m_dDelta = glass3::util::GlassMath::k_RadiansToDegrees
			* m_geoTTOrigin.delta(geo);

	return (T(m_dDelta));
}

// ---------------------------------------------------------T

double CTravelTime::T(double delta) {
	// bounds checks
	if(delta < m_dMinimumDistance || delta > m_dMaximumDistance) {
    return (k_dTravelTimeInvalid);
	}
	if(m_dDepth < m_dMinimumDepth || m_dDepth > m_dMaximumDepth) {
    return (k_dTravelTimeInvalid);
	}

	// Calculate travel time given delta in degrees
	// double depth = m_pDepthWarp->calculateGridPoint(m_dDepth);
	// double distance = m_pDistanceWarp->calculateGridPoint(delta);

	// compute travel time using bilinear interpolation
	double travelTime = bilinear(delta, m_dDepth);
	m_dDelta = delta;

	return (travelTime);
}

// ---------------------------------------------------------T
double CTravelTime::T(int deltaIndex, int depthIndex) {
	// bounds checks
	if ((deltaIndex < 0) || (deltaIndex >= m_iNumDistances)) {
		return (k_dTravelTimeInvalid);
	}
	if ((depthIndex < 0) || (depthIndex >= m_iNumDepths)) {
		return (k_dTravelTimeInvalid);
	}

	// get traveltime from travel time array
	double travelTime = m_pDepthTravelTimeArray[depthIndex * m_iNumDistances
			+ deltaIndex];

	return (travelTime);
}

// ---------------------------------------------------------Bilinear
double CTravelTime::bilinear(double distance, double depth) {
	double interpolationGrid[2][2];
	double travelTime;
	int startingDelta = static_cast<int>(distance);
	int startingDepth = static_cast<int>(depth);
	bool error = false;

	// generate interpolation grid
	for (int i = 0; i < 2; i++) {
		// calculate distance index
		int deltaIndex = startingDelta + i;

		for (int j = 0; j < 2; j++) {
			// calculate depth index
			int depthIndex = startingDepth + j;

			// get current travel time from travel time array
			double time = T(deltaIndex, depthIndex);

			// check current travel time
			if (time < 0.0) {
				error = true;
			}

			// store travel time in interpolation grid
			interpolationGrid[i][j] = time;
		}

		double s = distance - floor(distance);
		double t = depth - floor(depth);

		// compute overall travel time by interpolating grid
		travelTime = interpolationGrid[0][0] * (1.0f - s) * (1.0f - t)
				+ interpolationGrid[0][1] * (1.0f - s) * t
				+ interpolationGrid[1][0] * s * (1.0f - t)
				+ interpolationGrid[1][1] * s * t;
	}

	// check if we had errors
	if (error) {
		// no traveltime
		return (k_dTravelTimeInvalid);
	}

	return (travelTime);
}
}  // namespace traveltime
