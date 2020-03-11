#include <gtest/gtest.h>

int main(int argc, char **argv) {
	int rc = 0;
	try {
		::testing::InitGoogleTest(&argc, argv);
		rc = RUN_ALL_TESTS();
	} catch (...) {}

	return rc;
}