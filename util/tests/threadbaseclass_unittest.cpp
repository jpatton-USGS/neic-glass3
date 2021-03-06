#include <gtest/gtest.h>
#include <threadbaseclass.h>
#include <logger.h>
#include <string>

#define TESTTHREADNAME "threadbasestub"
#define TESTSLEEPTIME 50
#define TESTSLEEPTIMECHANGE 100
#define TESTSTOPCOUNT 10
#define WAITTIME 2

// stub class based on ThreadBaseClass for unit tests
class threadbasestub : public glass3::util::ThreadBaseClass {
 public:
	// basic constructor
	threadbasestub()
			: glass3::util::ThreadBaseClass() {
		runcount = 0;
		startcount = false;
		setHealthCheckInterval(1);
		kill = false;
	}

	// alternate constructor
	threadbasestub(std::string& name, int time)
			: glass3::util::ThreadBaseClass(name, time) {
		runcount = 0;
		startcount = false;
		setHealthCheckInterval(1);
		kill = false;
	}

	~threadbasestub() {
	}

	// flag to control when the work threads start counting
	bool startcount;

	// flag to shut down the work threads
	bool kill;

	// count the runs
	int runcount;

	// function to expose protected setThreadState() function from
	// ThreadBaseClass
	void testSetThreadState(glass3::util::ThreadState state) {
		setWorkThreadsState(state);
	}

 protected:
	// work function for tests
	glass3::util::WorkState work() override {
		// check for shutdown
		if (kill) {
			return (glass3::util::WorkState::Error);
		}

		// keep running until count
		// is 10
		if ((startcount) && (runcount < TESTSTOPCOUNT)) {
			// increase count
			runcount++;

			// work successful
			setThreadHealth();
			return (glass3::util::WorkState::OK);
		}

		// still alive
		setThreadHealth();
		return (glass3::util::WorkState::Idle);
	}
};

// tests to see if the threadbaseclass is functional
TEST(ThreadBaseClassTest, CombinedTest) {
	std::string name = std::string(TESTTHREADNAME);

	// create a threadbasestub
	threadbasestub * TestThreadBaseStub = new threadbasestub(name,
	TESTSLEEPTIME);

	// assert that class not started
	ASSERT_TRUE(TestThreadBaseStub->getWorkThreadsState() ==
			glass3::util::ThreadState::Initialized)<<
	"TestThreadBaseStub not started";

	// assert that healthCheck is true
	ASSERT_TRUE(TestThreadBaseStub->healthCheck())<<
	"TestThreadBaseStub healthCheck is true";

	// check sleep time
	ASSERT_EQ(TestThreadBaseStub->getSleepTime(), TESTSLEEPTIME)<<
	"Check sleep time";

	// change sleeptime
	TestThreadBaseStub->setSleepTime(TESTSLEEPTIMECHANGE);

	// check changed sleep time
	ASSERT_EQ(TestThreadBaseStub->getSleepTime(), TESTSLEEPTIMECHANGE)<<
	"Check changed sleep time";

	// check thread name
	std::string testthreadname = std::string(TESTTHREADNAME);
	ASSERT_STREQ(TestThreadBaseStub->getThreadName().c_str(),
			testthreadname.c_str())<< "check thread name";

	// start the thread
	ASSERT_TRUE(TestThreadBaseStub->start())<< "start was successful";

	// wait a little while
	std::this_thread::sleep_for(std::chrono::seconds(WAITTIME / 2));

	// assert that class running
	ASSERT_TRUE(TestThreadBaseStub->getWorkThreadsState() ==
			glass3::util::ThreadState::Started)<< "TestThreadBaseStub running";

	// assert that healthCheck is true
	ASSERT_TRUE(TestThreadBaseStub->healthCheck())<<
	"TestThreadBaseStub healthCheck is true";

	// tell stub to start counting
	TestThreadBaseStub->startcount = true;

	// wait a little while
	std::this_thread::sleep_for(std::chrono::seconds(WAITTIME));

	// assert that healthCheck is true
	ASSERT_TRUE(TestThreadBaseStub->healthCheck())<<
	"TestThreadBaseStub healthCheck is true";

	// check count
	ASSERT_EQ(TestThreadBaseStub->runcount, TESTSTOPCOUNT)<< "Check count";

	// stop the thread
	ASSERT_TRUE(TestThreadBaseStub->stop())<< "stop was successful";

	// wait a little while
	std::this_thread::sleep_for(std::chrono::seconds(WAITTIME / 2));

	// assert that class not started
	ASSERT_TRUE(TestThreadBaseStub->getWorkThreadsState() ==
			glass3::util::ThreadState::Stopped)<< "TestThreadBaseStub not started";

	// assert that healthCheck is false
	ASSERT_FALSE(TestThreadBaseStub->healthCheck())<<
	"TestThreadBaseStub healthCheck is false";

	// cleanup
	delete (TestThreadBaseStub);
}

// tests to see if the queue is functional
TEST(ThreadBaseClassTest, KillTest) {
	std::string name = std::string(TESTTHREADNAME);

	// create a threadbasestub
	threadbasestub * TestThreadBaseStub = new threadbasestub(name,
	TESTSLEEPTIME);

	// start the thread
	ASSERT_TRUE(TestThreadBaseStub->start())<< "start was successful";

	// wait a little while
	std::this_thread::sleep_for(std::chrono::seconds(WAITTIME));

	// assert that healthCheck is true
	ASSERT_TRUE(TestThreadBaseStub->healthCheck())<<
	"TestThreadBaseStub healthCheck is true";

	// kill it
	TestThreadBaseStub->kill = true;

	// wait a little while
	std::this_thread::sleep_for(std::chrono::seconds(WAITTIME));

	// assert that healthCheck is false
	ASSERT_FALSE(TestThreadBaseStub->healthCheck())<<
	"TestThreadBaseStub healthCheck is false";
}

// tests various failure conditions
TEST(ThreadBaseClassTest, FailTests) {
	// glass3::util::log_init("ThreadBaseClassTest", "debug", ".", true);

	// create a threadbasestub
	threadbasestub * TestThreadBaseStub = new threadbasestub();

	// assert we can't stop what we've not started
	ASSERT_FALSE(TestThreadBaseStub->stop())<<
	"TestThreadBaseStub stop is false";

	// start the thread
	ASSERT_TRUE(TestThreadBaseStub->start())<< "start was successful";
	ASSERT_FALSE(TestThreadBaseStub->start())<<
	"second start was not successful";

	// allocation test
	TestThreadBaseStub->testSetThreadState(glass3::util::ThreadState::Stopped);
	ASSERT_FALSE(TestThreadBaseStub->start())<<
	"second start (allocation) was not successful";

	// kill it
	TestThreadBaseStub->kill = true;

	// wait a little while
	std::this_thread::sleep_for(std::chrono::seconds(WAITTIME));
	ASSERT_FALSE(TestThreadBaseStub->stop())<<
	"stop after fail was not successful";

	TestThreadBaseStub->testSetThreadState(glass3::util::ThreadState::Started);
	TestThreadBaseStub->setThreadHealth(false);
	ASSERT_FALSE(TestThreadBaseStub->healthCheck())<<
	"TestThreadBaseStub healthCheck is false";
}
