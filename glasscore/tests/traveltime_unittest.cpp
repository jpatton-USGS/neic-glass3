#include <gtest/gtest.h>

#include <string>

#include <logger.h>

#include "TravelTime.h"

#define TESTPATH "testdata"
#define PHASE "P"
#define PHASEFILENAME "P.ttt"

#define NDISTANCE 550
#define NDEPTH 105
#define NTIME 105

#define MINDIST 0
#define MAXDIST 180

#define MINDEPTH 0
#define MAXDEPTH 900

#define MINTIME 0
#define MAXTIME 180

#define DELTADIST 14
#define DELTADEPTH 15
#define DELTATIME 16

#define LATITUDE 0.0
#define LONGITUDE 0.0
#define DEPTH 50.0
#define DISTANCE 50.0
#define TIME 529.2
#define TIME2 50.553
#define BILINEAR 50.553

// tests to see if the traveltime can be constructed
TEST(TravelTimeTest, Construction) {
	glass3::util::Logger::disable();

	// construct a traveltime
	traveltime::CTravelTime traveltime;

	// m_iNumDistances
	ASSERT_EQ(0, traveltime.m_iNumDistances)<< "m_iNumDistances Check";

	// m_dMinimumDistance
	ASSERT_EQ(0, traveltime.m_dMinimumDistance)<< "m_dMinimumDistance Check";

	// m_dMaximumDistance
	ASSERT_EQ(0, traveltime.m_dMaximumDistance)<< "m_dMaximumDistance Check";

	// m_iNumDepths
	ASSERT_EQ(0, traveltime.m_iNumDepths)<< "m_iNumDepths Check";

	// m_dMinimumDepth
	ASSERT_EQ(0, traveltime.m_dMinimumDepth)<< "m_dMinimumDepth Check";

	// m_dMaximumDepth
	ASSERT_EQ(0, traveltime.m_dMaximumDepth)<< "m_dMaximumDepth Check";

	// m_iNumTimes
	ASSERT_EQ(0, traveltime.m_iNumTimes)<< "m_iNumTimes Check";

	// m_dMinimumTime
	ASSERT_EQ(0, traveltime.m_dMinimumTime)<< "m_dMinimumTime Check";

	// m_dMaximumTime
	ASSERT_EQ(0, traveltime.m_dMaximumTime)<< "m_dMaximumTime Check";

	// m_dDistanceDelta
	ASSERT_EQ(0, traveltime.m_dDistanceDelta)<< "m_dDistanceDelta Check";

	// m_dDepthDelta
	ASSERT_EQ(0, traveltime.m_dDepthDelta)<< "m_dDepthDelta Check";

	// m_dTimeDelta
	ASSERT_EQ(0, traveltime.m_dTimeDelta)<< "m_dTimeDelta Check";

	// dDepth
	ASSERT_EQ(0, traveltime.m_dDepth)<< "Depth Check";

	// dDelta
	ASSERT_EQ(0, traveltime.m_dDelta)<< "Delta Check";

	// pointers
	ASSERT_EQ(NULL, traveltime. m_pDepthTravelTimeArray)<< "m_pDepthTravelTimeArray " // NOLINT
		"null";
	ASSERT_EQ(NULL, traveltime.m_pDepthDistanceArray)<< "pDepthDistanceArray null";
}

// tests to see if the traveltime can be setup
TEST(TravelTimeTest, Setup) {
	glass3::util::Logger::disable();

	std::string phasefile = "./" + std::string(TESTPATH) + "/"
			+ std::string(PHASEFILENAME);
	std::string phasename = std::string(PHASE);

	// construct a traveltime
	traveltime::CTravelTime traveltime;

	// setup
	traveltime.setup(phasename, phasefile);

	// phase name
	ASSERT_STREQ(traveltime.m_sPhase.c_str(), phasename.c_str());

	// m_iNumDistances
	ASSERT_EQ(NDISTANCE, traveltime.m_iNumDistances)<< "m_iNumDistances Check";

	// m_dMinimumDistance
	ASSERT_EQ(MINDIST, traveltime.m_dMinimumDistance)<< "m_dMinimumDistance Check";

	// m_dMaximumDistance
	ASSERT_EQ(MAXDIST, traveltime.m_dMaximumDistance)<< "m_dMaximumDistance Check";

	// m_iNumDepths
	ASSERT_EQ(NDEPTH, traveltime.m_iNumDepths)<< "m_iNumDepths Check";

	// m_dMinimumDepth
	ASSERT_EQ(MINDEPTH, traveltime.m_dMinimumDepth)<< "m_dMinimumDepth Check";

	// m_dMaximumDepth
	ASSERT_EQ(MAXDEPTH, traveltime.m_dMaximumDepth)<< "m_dMaximumDepth Check";

	// m_iNumTimes
	ASSERT_EQ(NTIME, traveltime.m_iNumTimes)<< "m_iNumTimes Check";

	// m_dMinimumTime
	ASSERT_EQ(MINTIME, traveltime.m_dMinimumTime)<< "m_dMinimumTime Check";

	// m_dMaximumTime
	ASSERT_EQ(MAXTIME, traveltime.m_dMaximumTime)<< "m_dMaximumTime Check";

	// m_dDistanceDelta
	ASSERT_EQ(DELTADIST, traveltime.m_dDistanceDelta)<< "m_dDistanceDelta Check";

	// m_dDepthDelta
	ASSERT_EQ(DELTADEPTH, traveltime.m_dDepthDelta)<< "m_dDepthDelta Check";

	// m_dTimeDelta
	ASSERT_EQ(DELTATIME, traveltime.m_dTimeDelta)<< "m_dTimeDelta Check";

	// dDepth
	ASSERT_EQ(DEPTH, traveltime.m_dDepth)<< "Depth Check";

	// dDelta
	ASSERT_EQ(0, traveltime.m_dDelta)<< "Delta Check";

	// pointers
	ASSERT_TRUE(NULL != traveltime. m_pDepthTravelTimeArray)<< "m_pDepthTravelTimeArray not " // NOLINT
			"null";
	ASSERT_TRUE(NULL != traveltime.m_pDepthDistanceArray)<< "pDepthDistanceArray "
			"not null";
}

// tests the time warp copy constructor
TEST(TravelTimeTest, Copy) {
	glass3::util::Logger::disable();

	std::string phasefile = "./" + std::string(TESTPATH) + "/"
			+ std::string(PHASEFILENAME);
	std::string phasename = std::string(PHASE);

	// construct a second traveltime
	traveltime::CTravelTime traveltime2;

	// setup
	traveltime2.setup(phasename, phasefile);

	// set origin
	traveltime2.setTTOrigin(LATITUDE, LONGITUDE, DEPTH);

	// copy
	traveltime::CTravelTime traveltime1(traveltime2);

	// phase name
	ASSERT_STREQ(traveltime1.m_sPhase.c_str(), phasename.c_str());

	// nDistanceWarp
	ASSERT_EQ(NDISTANCEWARP, traveltime1.m_iNumDistanceWarp)<< "nDistanceWarp Check";

	// nDepthWarp
	ASSERT_EQ(NDEPTHWARP, traveltime1.m_iNumDepthWarp)<< "nDepthWarp Check";

	// dDepth
	ASSERT_NEAR(DEPTH, traveltime1.m_dDepth, 0.001)<< "Depth Check";

	// dDelta
	ASSERT_EQ(DISTANCE, traveltime1.m_dDelta)<< "Delta Check";

	// pointers
	ASSERT_TRUE(NULL != traveltime1.m_pDistanceWarp)<< "pDistanceWarp not null";
	ASSERT_TRUE(NULL != traveltime1.m_pDepthWarp)<< "pDepthWarp not null";
	ASSERT_TRUE(NULL != traveltime1. m_pTravelTimeArray)<< "pTravelTimeArray not "
			"null";
	ASSERT_TRUE(NULL != traveltime1.m_pDepthDistanceArray)<< "pDepthDistanceArray "
			"not null";
	ASSERT_TRUE(NULL != traveltime1.m_pPhaseArray)<< "pPhaseArray not null";
}

// tests traveltime operations
TEST(TravelTimeTest, Operations) {
	glass3::util::Logger::disable();

	std::string phasefile = "./" + std::string(TESTPATH) + "/"
			+ std::string(PHASEFILENAME);
	std::string phasename = std::string(PHASE);

	// construct a traveltime
	traveltime::CTravelTime traveltime;

	// setup
	traveltime.setup(phasename, phasefile);

	// set origin
	traveltime.setTTOrigin(LATITUDE, LONGITUDE, DEPTH);

	// T(delta)
	ASSERT_NEAR(TIME, traveltime.T(DISTANCE), 0.001)<< "T(delta) Check";

	// dDepth
	ASSERT_NEAR(DEPTH, traveltime.m_dDepth, 0.001)<< "Depth Check";

	// dDelta
	ASSERT_NEAR(DISTANCE, traveltime.m_dDelta, 0.001)<< "Delta Check";

	glass3::util::Geo testGeo;
	testGeo.setGeographic(LATITUDE, LONGITUDE + DISTANCE, DEPTH);

	// T(geo)
	ASSERT_NEAR(TIME, traveltime.T(&testGeo), 0.001)<< "T(geo) Check";

	// dDepth
	ASSERT_NEAR(DEPTH, traveltime.m_dDepth, 0.001)<< "Depth Check";

	// dDelta
	ASSERT_NEAR(DISTANCE, traveltime.m_dDelta, 0.001)<< "Delta Check";

	// T(delta, distance)
	ASSERT_NEAR(TIME2, traveltime.T(DISTANCE,DEPTH), 0.001)<< "T(delta, distance) Check"; // NOLINT

	// bilinear
	ASSERT_NEAR(BILINEAR, traveltime.bilinear(DISTANCE,DEPTH), 0.001)<< "bilinear Check"; // NOLINT
}
