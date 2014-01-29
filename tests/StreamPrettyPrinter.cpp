#include <fstream>
#include <iostream>
#include <string>

#include "gtest/gtest.h"
#include "gtest/internal/gtest-internal.h"

#include "tests/StreamPrettyPrinter.hpp"

using namespace std;
using namespace testing;
using namespace testing::internal;

// copied from gtest.cc
// Formats a countable noun.  Depending on its quantity, either the
// singular form or the plural form is used. e.g.
//
// FormatCountableNoun(1, "formula", "formuli") returns "1 formula".
// FormatCountableNoun(5, "book", "books") returns "5 books".
static std::string FormatCountableNoun(int count,
        const char * singular_form,
        const char * plural_form) {
    return internal::StreamableToString(count) + " " +
        (count == 1 ? singular_form : plural_form);
}

// Formats the count of tests.
static std::string FormatTestCount(int test_count) {
    return FormatCountableNoun(test_count, "test", "tests");
}

// Formats the count of test cases.
static std::string FormatTestCaseCount(int test_case_count) {
    return FormatCountableNoun(test_case_count, "test case", "test cases");
}

// Converts a TestPartResult::Type enum to human-friendly string
// representation.  Both kNonFatalFailure and kFatalFailure are translated
// to "Failure", as the user usually doesn't care about the difference
// between the two when viewing the test result.
static const char * TestPartResultTypeToString(TestPartResult::Type type) {
    switch (type) {
        case TestPartResult::kSuccess:
            return "Success";

        case TestPartResult::kNonFatalFailure:
        case TestPartResult::kFatalFailure:
#ifdef _MSC_VER
            return "error: ";
#else
            return "Failure\n";
#endif
        default:
            return "Unknown result type";
    }
}

static void PrintFullTestCommentIfPresent(
        ostream &out, const TestInfo& test_info) {
    const char* const type_param = test_info.type_param();
    const char* const value_param = test_info.value_param();

    if (type_param != NULL || value_param != NULL) {
        out << ", where ";
        if (type_param != NULL) {
            out << "TypeParam = " << type_param;
            if (value_param != NULL)
                out << " and ";
        }
        if (value_param != NULL) {
            out << "GetParam() = " << value_param;
        }
    }
}

// Prints a TestPartResult to an std::string.
static std::string PrintTestPartResultToString(
        const TestPartResult& test_part_result) {
    return (Message()
            << internal::FormatFileLocation(test_part_result.file_name(),
                test_part_result.line_number())
            << " " << TestPartResultTypeToString(test_part_result.type())
            << test_part_result.message()).GetString();
}

// Prints a TestPartResult.
static void PrintTestPartResult(
        ostream &out, const TestPartResult& test_part_result) {
    const std::string& result =
        PrintTestPartResultToString(test_part_result);
    out << result << endl;
    // If the test program runs in Visual Studio or a debugger, the
    // following statements add the test part result message to the Output
    // window such that the user can double-click on it to jump to the
    // corresponding source code location; otherwise they do nothing.
#if GTEST_OS_WINDOWS && !GTEST_OS_WINDOWS_MOBILE
    // We don't call OutputDebugString*() on Windows Mobile, as printing
    // to stdout is done by OutputDebugString() there already - we don't
    // want the same message printed twice.
    ::OutputDebugStringA(result.c_str());
    ::OutputDebugStringA("\n");
#endif
}


// Fired before each iteration of tests starts.
void StreamPrettyPrinter::OnTestIterationStart(
        const UnitTest& unit_test, int iteration) {
    if (GTEST_FLAG(repeat) != 1)
        out << endl
            << "Repeating all tests (iteration " << iteration + 1 << ") . . ."
            << endl
            << endl;

    const char* const filter = GTEST_FLAG(filter).c_str();

    // Prints the filter if it's not *.  This reminds the user that some
    // tests may be skipped.
    if (!String::CStringEquals(filter, "*")) {
        out << "Note: " << GTEST_NAME_ << " filter = " << filter << endl;
    }

    if (GTEST_FLAG(shuffle)) {
        out << "Note: Randomizing tests' orders with a seed of "
            << unit_test.random_seed()
            << " ."
            << endl;
    }

    out << "[==========] ";
    out << "Running "
        <<  FormatTestCount(unit_test.test_to_run_count())
        << " from "
        <<  FormatTestCaseCount(unit_test.test_case_to_run_count())
        << "."
        << endl;
}

void StreamPrettyPrinter::OnEnvironmentsSetUpStart(
        const UnitTest& /*unit_test*/) {
    out << "[----------] ";
    out << "Global test environment set-up." << endl;
}

void StreamPrettyPrinter::OnTestCaseStart(const TestCase& test_case) {
    const std::string counts =
        FormatCountableNoun(test_case.test_to_run_count(), "test", "tests");
    out << "[----------] ";
    out << counts << " from " << test_case.name();
    if (test_case.type_param() == NULL) {
        out << endl;
    } else {
        out << ", where TypeParam = " << test_case.type_param() << endl;
    }
}

void StreamPrettyPrinter::OnTestStart(const TestInfo& test_info) {
    out << "[ RUN      ] ";
    PrintTestName(test_info.test_case_name(), test_info.name());
    out << endl;
}

// Called after an assertion failure.
void StreamPrettyPrinter::OnTestPartResult(
        const TestPartResult& result) {
    // If the test part succeeded, we don't need to do anything.
    if (result.type() == TestPartResult::kSuccess)
        return;

    // Print failure message from the assertion (e.g. expected this and got that).
    PrintTestPartResult(out, result);
    fflush(stdout);
}

void StreamPrettyPrinter::OnTestEnd(const TestInfo& test_info) {
    if (test_info.result()->Passed()) {
        out << "[       OK ] ";
    } else {
        out << "[  FAILED  ] ";
    }
    PrintTestName(test_info.test_case_name(), test_info.name());
    if (test_info.result()->Failed())
        PrintFullTestCommentIfPresent(out, test_info);

    if (GTEST_FLAG(print_time)) {
        out << " ("
            << test_info.result()->elapsed_time()
            << " ms)" << endl;
    } else {
        out << endl;
    }
}

void StreamPrettyPrinter::OnTestCaseEnd(const TestCase& test_case) {
    if (!GTEST_FLAG(print_time)) return;

    const std::string counts =
        FormatCountableNoun(test_case.test_to_run_count(), "test", "tests");
    out << "[----------] "
        << counts
        << " from "
        << test_case.name()
        << " ("
        << internal::StreamableToString(test_case.elapsed_time())
        << " ms total)"
        << endl
        << endl;
}

void StreamPrettyPrinter::OnEnvironmentsTearDownStart(
        const UnitTest& /*unit_test*/) {
    out << "[----------] "
        << "Global test environment tear-down"
        << endl;
}

// Internal helper for printing the list of failed tests.
void StreamPrettyPrinter::PrintFailedTests(const UnitTest& unit_test) {
    const int failed_test_count = unit_test.failed_test_count();
    if (failed_test_count == 0) {
        return;
    }

    for (int i = 0; i < unit_test.total_test_case_count(); ++i) {
        const TestCase& test_case = *unit_test.GetTestCase(i);
        if (!test_case.should_run() || (test_case.failed_test_count() == 0)) {
            continue;
        }
        for (int j = 0; j < test_case.total_test_count(); ++j) {
            const TestInfo& test_info = *test_case.GetTestInfo(j);
            if (!test_info.should_run() || test_info.result()->Passed()) {
                continue;
            }
            out << "[  FAILED  ] ";
            out << test_case.name() << '.' << test_info.name();
            PrintFullTestCommentIfPresent(out, test_info);
            out << endl;
        }
    }
}

void StreamPrettyPrinter::OnTestIterationEnd(const UnitTest& unit_test,
        int /*iteration*/) {
    out << "[==========] ";
    out <<  FormatTestCount(unit_test.test_to_run_count())
        << " from "
        <<  FormatTestCaseCount(unit_test.test_case_to_run_count())
        << " ran.";
    if (GTEST_FLAG(print_time)) {
        out << " (" << unit_test.elapsed_time() << " ms total)";
    }
    out << endl;
    out << "[  PASSED  ] ";
    out << FormatTestCount(unit_test.successful_test_count()) << endl;

    int num_failures = unit_test.failed_test_count();
    if (!unit_test.Passed()) {
        const int failed_test_count = unit_test.failed_test_count();
        out << "[  FAILED  ] ";
        out << FormatTestCount(failed_test_count) << ", listed below:" << endl;
        PrintFailedTests(unit_test);
        out << endl
            << num_failures
            << " FAILED "
            << (num_failures == 1 ? "TEST" : "TESTS")
            << endl;
    }

    int num_disabled = unit_test.reportable_disabled_test_count();
    if (num_disabled && !GTEST_FLAG(also_run_disabled_tests)) {
        if (!num_failures) {
            out << endl;  // Add a spacer if no FAILURE banner is displayed.
        }
        out << "  YOU HAVE " << num_disabled << " DISABLED "
            << (num_disabled == 1 ? "TEST" : "TESTS")
            << endl
            << endl;
    }
    // Ensure that Google Test output is printed before, e.g., heapchecker output.
    fflush(stdout);
}

// End StreamPrettyPrinter

