#include "TTT.h"
#include <logger.h>
#include <geo.h>
#include <fstream>
#include <string>
#include <cmath>
#include "TravelTime.h"

namespace traveltime {

// constants
constexpr double CTTT::k_dTTTooLargeToBeValid;
const int CTTT::k_iMaximumNumberOfTravelTimes;

// ---------------------------------------------------------CTTT
CTTT::CTTT() {
	clear();
}

// ---------------------------------------------------------CTTT
CTTT::CTTT(const CTTT &ttt) {
	clear();

	m_iNumTravelTimes = ttt.m_iNumTravelTimes;
	m_geoTTOrigin = ttt.m_geoTTOrigin;

	for (int i = 0; i < ttt.m_iNumTravelTimes; i++) {
		if (ttt.m_pTravelTimes[i] != NULL) {
			m_pTravelTimes[i] = new CTravelTime(*ttt.m_pTravelTimes[i]);
		} else {
			m_pTravelTimes[i] = NULL;
		}

		m_adMinimumAssociationValues[i] = ttt.m_adMinimumAssociationValues[i];
		m_adMaximumAssociationValues[i] = ttt.m_adMaximumAssociationValues[i];
	}
}

// ---------------------------------------------------------~CTTT
CTTT::~CTTT() {
	for (int i = 0; i < m_iNumTravelTimes; i++) {
		delete (m_pTravelTimes[i]);
	}
}

void CTTT::clear() {
	m_iNumTravelTimes = 0;
	m_geoTTOrigin.clear();

	for (int i = 0; i < k_iMaximumNumberOfTravelTimes; i++) {
		m_pTravelTimes[i] = NULL;
		m_adMinimumAssociationValues[i] = CTravelTime::k_dTravelTimeInvalid;
		m_adMaximumAssociationValues[i] = CTravelTime::k_dTravelTimeInvalid;
	}
}

// -------------------------------------------------------writeToFiles
void CTTT::writeToFiles(std::string outDir, double depth) {
	// write out each ttfile we know about as csv
	for (int i = 0; i < k_iMaximumNumberOfTravelTimes; i++) {
		if (m_pTravelTimes[i] != NULL) {
			std::string fileName = outDir + "/" + m_pTravelTimes[i]->m_sPhase
				+ "_" + std::to_string(depth) + ".csv";
			m_pTravelTimes[i]->writeToFile(fileName, depth);
		}
	}
}

// ---------------------------------------------------------addPhase
// Add phase to list to be calculated
bool CTTT::addPhase(std::string phase, double *assocRange,
					std::string file, bool useForLocation,
					double minPublishable, double maxPublishable) {
	// bounds check
	if ((m_iNumTravelTimes + 1) > k_iMaximumNumberOfTravelTimes) {
		glass3::util::Logger::log(
				"error",
				"CTTT::addPhase: Maximum number of phases ("
						+ std::to_string(k_iMaximumNumberOfTravelTimes) + ") reached");
		return (false);
	}

	// create and setup traveltime from phase
	CTravelTime *trv = new CTravelTime(useForLocation, minPublishable,
									   maxPublishable);
	trv->setup(phase, file);

	// add traveltime to list
	m_pTravelTimes[m_iNumTravelTimes] = trv;

	// set up association range
	if (assocRange != NULL) {
		m_adMinimumAssociationValues[m_iNumTravelTimes] = assocRange[0];
		m_adMaximumAssociationValues[m_iNumTravelTimes] = assocRange[1];
	}

	m_iNumTravelTimes++;

	return (true);
}

// ---------------------------------------------------------setTTOrigin
void CTTT::setTTOrigin(double lat, double lon, double z) {
	// this should go ahead and update the CGeo
	m_geoTTOrigin.setGeographic(lat, lon,
								glass3::util::Geo::k_EarthRadiusKm - z);
}

// ---------------------------------------------------------setTTOrigin
void CTTT::setTTOrigin(const glass3::util::Geo &geoOrigin) {
	m_geoTTOrigin = geoOrigin;
}

// ---------------------------------------------------------T
double CTTT::T(glass3::util::Geo *geo, std::string phase) {
	// Calculate travel time from distance in degrees
	// for each phase
	for (int i = 0; i < m_iNumTravelTimes; i++) {
		// is this the phase we're looking for
		if (m_pTravelTimes[i]->m_sPhase == phase) {
			// set origin
			m_pTravelTimes[i]->setTTOrigin(m_geoTTOrigin);

			// get travel time and phase
			double traveltime = m_pTravelTimes[i]->T(geo);
			m_sPhase = phase;
			m_bUseForLocations = m_pTravelTimes[i]->m_bUseForLocations;

			if ((m_pTravelTimes[i]->m_dDelta >=
					m_pTravelTimes[i]->m_dMinDeltaPublishable) &&
				(m_pTravelTimes[i]->m_dDelta <=
					m_pTravelTimes[i]->m_dMaxDeltaPublishable)) {
				m_bPublishable = true;
			} else {
				m_bPublishable = false;
			}

			return (traveltime);
		}
	}

	// no valid travel time
	m_sPhase = "?";
	m_bUseForLocations = false;
	m_bPublishable = false;
	return (CTravelTime::k_dTravelTimeInvalid);
}

// ---------------------------------------------------------T
double CTTT::Td(double delta, std::string phase, double depth) {
	m_geoTTOrigin.m_dGeocentricRadius = glass3::util::Geo::k_EarthRadiusKm
			- depth;
	// Calculate time from delta (degrees) and depth
	// for each phase
	for (int i = 0; i < m_iNumTravelTimes; i++) {
		// is this the phase we're looking for
		if (m_pTravelTimes[i]->m_sPhase == phase) {
			// set origin and depth
			m_pTravelTimes[i]->setTTOrigin(m_geoTTOrigin);

			// get travel time and phase
			double traveltime = m_pTravelTimes[i]->T(delta);
			m_sPhase = phase;
			m_bUseForLocations = m_pTravelTimes[i]->m_bUseForLocations;

			if ((m_pTravelTimes[i]->m_dDelta >=
					m_pTravelTimes[i]->m_dMinDeltaPublishable) &&
				(m_pTravelTimes[i]->m_dDelta <=
					m_pTravelTimes[i]->m_dMaxDeltaPublishable)) {
				m_bPublishable = true;
			} else {
				m_bPublishable = false;
			}

			return (traveltime);
		}
	}

	// no valid travel time
	m_sPhase = "?";
	m_bUseForLocations = false;
	m_bPublishable = false;
	return (CTravelTime::k_dTravelTimeInvalid);
}

// ---------------------------------------------------------T
double CTTT::T(double delta, std::string phase) {
	// Calculate time from delta (degrees)
	// for each phase
	for (int i = 0; i < m_iNumTravelTimes; i++) {
		// is this the phase we're looking for
		if (m_pTravelTimes[i]->m_sPhase == phase) {
			// set origin
			m_pTravelTimes[i]->setTTOrigin(m_geoTTOrigin);

			// get travel time and phase
			double traveltime = m_pTravelTimes[i]->T(delta);
			m_sPhase = phase;
			m_bUseForLocations = m_pTravelTimes[i]->m_bUseForLocations;

			if ((m_pTravelTimes[i]->m_dDelta >=
					m_pTravelTimes[i]->m_dMinDeltaPublishable) &&
				(m_pTravelTimes[i]->m_dDelta <=
					m_pTravelTimes[i]->m_dMaxDeltaPublishable)) {
				m_bPublishable = true;
			} else {
				m_bPublishable = false;
			}

			return (traveltime);
		}
	}

	// no valid travel time
	m_sPhase = "?";
	m_bUseForLocations = false;
	m_bPublishable = false;
	return (CTravelTime::k_dTravelTimeInvalid);
}

// ---------------------------------------------------------testTravelTimes
/* double CTTT::testTravelTimes(std::string phase) {
	// Calculate time from delta (degrees)
	double time = 0;
	std::ofstream outfile;
	std::string filename = phase + "_travel_time_Z_0.txt";
	outfile.open(filename, std::ios::out);

	for (double i = 0; i < 5.; i += 0.005) {
		time = Td(i, phase, 0.0);
		outfile << std::to_string(i) << ", " << std::to_string(time) << "\n";
	}
	outfile.close();

	filename = phase + "_travel_time_Z_50.txt";
	outfile.open(filename, std::ios::out);

	for (double i = 0; i < 5.; i += 0.005) {
		time = Td(i, phase, 50.);
		outfile << std::to_string(i) << ", " << std::to_string(time) << "\n";
	}
	outfile.close();

	filename = phase + "_travel_time_Z_100.txt";
	outfile.open(filename, std::ios::out);

	for (double i = 0; i < 5.; i += 0.005) {
		time = Td(i, phase, 100.);
		outfile << std::to_string(i) << ", " << std::to_string(time) << "\n";
	}
	outfile.close();

	return (1.0);
} */

// ---------------------------------------------------------T
double CTTT::T(glass3::util::Geo *geo, double tObserved) {
	// Find Phase with least residual, returns time

	double bestTraveltime = CTravelTime::k_dTravelTimeInvalid;
	std::string bestPhase = "?";
	bool useForLocations = true;
	bool publishable = true;
	double bestResidual = k_dTTTooLargeToBeValid;

	// for each phase
	for (int i = 0; i < m_iNumTravelTimes; i++) {
		// get current aTrv
		CTravelTime * aTrv = m_pTravelTimes[i];

		// set origin
		aTrv->setTTOrigin(m_geoTTOrigin);

		// get traveltime
		double traveltime = aTrv->T(geo);

		// check traveltime
		if (traveltime <= CTravelTime::k_dTravelTimeInvalid) {
			continue;
		}

		// check to see if phase is associable
		// based on minimum assoc distance, if present
		if (m_adMinimumAssociationValues[i] >= 0) {
			if (aTrv->m_dDelta < m_adMinimumAssociationValues[i]) {
				// this phase is not associable  at this distance
				continue;
			}
		}

		// check to see if phase is associable
		// based on maximum assoc distance, if present
		if (m_adMaximumAssociationValues[i] >= 0) {
			if (aTrv->m_dDelta > m_adMaximumAssociationValues[i]) {
				// this phase is not associable  at this distance
				continue;
			}
		}

		// compute residual
		double residual = std::abs(tObserved - traveltime);

		// check to see if this residual is better than the previous
		// best
		if (residual < bestResidual) {
			// this is the new best travel time
			bestResidual = residual;
			bestPhase = aTrv->m_sPhase;
			bestTraveltime = traveltime;
			useForLocations = aTrv->m_bUseForLocations;

			if ((aTrv->m_dDelta >= aTrv->m_dMinDeltaPublishable) &&
				(aTrv->m_dDelta <= aTrv->m_dMaxDeltaPublishable)) {
				publishable = true;
			} else {
				publishable = false;
			}
		}
	}

	// check to see if minimum residual is valid
	if (bestResidual < k_dTTTooLargeToBeValid) {
		m_sPhase = bestPhase;
		m_bUseForLocations = useForLocations;
		m_bPublishable = publishable;

		return (bestTraveltime);
	}

	// no valid travel time
	m_sPhase = "?";
	m_bUseForLocations = false;
	m_bPublishable = false;
	return (CTravelTime::k_dTravelTimeInvalid);
}
}  // namespace traveltime
