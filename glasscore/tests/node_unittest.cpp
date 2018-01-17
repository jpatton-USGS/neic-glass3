#include <gtest/gtest.h>
#include <memory>
#include <string>
#include "Node.h"
#include "Site.h"
#include "Logit.h"

// test data
#define NAME "testNode"
#define LATITUDE 45.90
#define LONGITUDE -112.79
#define DEPTH 10
#define RESOLUTION 25
#define NODEID "1234567890"

#define SITEJSON "{\"Cmd\":\"Site\",\"Elv\":2326.000000,\"Lat\":45.822170,\"Lon\":-112.451000,\"Site\":\"LRM.EHZ.MB.--\",\"Use\":true}"  // NOLINT
#define TRAVELTIME 122

// NOTE: Need to consider testing nucleate, but that would need a much more
// involved set of real data, and possibly a glass refactor to better support
// unit testing

// tests to see if the node can be constructed
TEST(NodeTest, Construction) {
	glassutil::CLogit::disable();

	// construct a node
	glasscore::CNode * testNode = new glasscore::CNode(std::string(NAME),
	LATITUDE,
														LONGITUDE, DEPTH,
														RESOLUTION,
														std::string(NODEID));

	// name
	ASSERT_STREQ(std::string(NAME).c_str(), testNode->sName.c_str())<<
			"Node Name Matches";

	// latitude
	ASSERT_EQ(LATITUDE, testNode->dLat)<< "Node Latitude Check";

	// longitude
	ASSERT_EQ(LONGITUDE, testNode->dLon)<< "Node Longitude Check";

	// depth
	ASSERT_EQ(DEPTH, testNode->dZ)<< "Node Depth Check";

	// resolution
	ASSERT_EQ(RESOLUTION, testNode->dResolution)<< "Node Resolution Check";

	// id
	ASSERT_STREQ(std::string(NODEID).c_str(), testNode->sPid.c_str())<<
			"Node ID Matches";

	// site list
	int expectedSize = 0;
	ASSERT_EQ(expectedSize, testNode->getSiteLinksCount())<< "sitelist empty";
}

// tests to see if sites can be added to the node
TEST(NodeTest, SiteOperations) {
	glassutil::CLogit::disable();

	// construct a node
	glasscore::CNode * testNode = new glasscore::CNode(std::string(NAME),
	LATITUDE,
														LONGITUDE, DEPTH,
														RESOLUTION,
														std::string(NODEID));

	// create json object from the string
	json::Object siteJSON = json::Deserialize(std::string(SITEJSON));

	// construct sites using JSON objects
	glasscore::CSite * testSite = new glasscore::CSite(&siteJSON, NULL);

	// create new shared pointer to the node
	std::shared_ptr<glasscore::CNode> sharedTestNode(testNode);

	// create new shared pointer to the site
	std::shared_ptr<glasscore::CSite> sharedTestSite(testSite);

	// add the site to the node
	testNode->linkSite(sharedTestSite, sharedTestNode, TRAVELTIME);

	// check to see if the site was added
	int expectedSize = 1;
	ASSERT_EQ(expectedSize, testNode->getSiteLinksCount())<< "node has one site";

	// test clearing the node site links
	testNode->clearSiteLinks();

	// check to see if the site list was cleared
	expectedSize = 0;
	ASSERT_EQ(expectedSize, testNode->getSiteLinksCount())<< "sitelist cleared";
}