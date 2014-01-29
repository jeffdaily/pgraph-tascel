#ifndef _PGRAPH_STREAMPRETTYPRINTER_H_
#define _PGRAPH_STREAMPRETTYPRINTER_H_

#include <fstream>
#include <iostream>
#include <string>

#include "gtest/gtest.h"

using namespace std;
using namespace testing;

// This class implements the TestEventListener interface.
//
// Class StreamPrettyPrinter is copyable.
class StreamPrettyPrinter : public TestEventListener {
    public:
        StreamPrettyPrinter(ostream &out)
            :   TestEventListener()
            ,   out(out)
        {}

        virtual ~StreamPrettyPrinter() {}

        // The following methods override what's in the TestEventListener class.
        virtual void OnTestProgramStart(const UnitTest& /*unit_test*/) {}
        virtual void OnTestIterationStart(const UnitTest& unit_test, int iteration);
        virtual void OnEnvironmentsSetUpStart(const UnitTest& unit_test);
        virtual void OnEnvironmentsSetUpEnd(const UnitTest& /*unit_test*/) {}
        virtual void OnTestCaseStart(const TestCase& test_case);
        virtual void OnTestStart(const TestInfo& test_info);
        virtual void OnTestPartResult(const TestPartResult& result);
        virtual void OnTestEnd(const TestInfo& test_info);
        virtual void OnTestCaseEnd(const TestCase& test_case);
        virtual void OnEnvironmentsTearDownStart(const UnitTest& unit_test);
        virtual void OnEnvironmentsTearDownEnd(const UnitTest& /*unit_test*/) {}
        virtual void OnTestIterationEnd(const UnitTest& unit_test, int iteration);
        virtual void OnTestProgramEnd(const UnitTest& /*unit_test*/) {}

    private:
        void PrintTestName(const string &test_case, const string &test) {
            out << test_case << '.' << test;
        }

        void PrintFailedTests(const UnitTest& unit_test);

        ostream &out;
};

#endif /* _PGRAPH_STREAMPRETTYPRINTER_H_ */
