/*****************************************
 * This file is documented for Doxygen.
 * If you modify this file please update
 * the comments so that Doxygen will still
 * be able to work.
 ****************************************/
/**
 * \file
 * \brief stringutil.h
 *
 * stringutil.h is a set of functions that manage string functions not provided
 * by the standard template library, including splitting strings and replacing
 * characters within strings.
 */
#ifndef STRINGUTIL_H
#define STRINGUTIL_H

#include <string>
#include <vector>

namespace glass3 {
namespace util {
/**
 * \brief split a string into substrings
 *
 * Split a string into substrings using the provided delimiter,
 * creating a std::vector to hold the results
 *
 * \param sInput - A std::string containing the string to split
 * \param cDelimiter - A char containing the the delimiter to split with.
 * \return returns a std::vector containing the split std::string elements
 */
std::vector<std::string> split(const std::string &sInput, char cDelimiter);

/**
 * \brief remove all instances of characters from a string
 *
 * Remove all instances of the set of characters stored in chars from the given
 * string.
 * This is NOT removing a substring from a string.
 *
 * \param sInput - A std::string containing the string to remove characters from
 * \param sRemoveChars - A std::string containing the individual characters to
 * remove.
 * \return returns a std::string containing modified string.
 */
std::string& removeChars(std::string& sInput, const std::string& sRemoveChars);  // NOLINT

/**
 * \brief convert provided double to a string with desired precision
 *
 * Convert the provided double value into a string while truncateing it to the
 * provided precision
 *
 * \param doubleVal - A double containing the value to convert
 * \param precisionVal - An int containing the desired precision to truncate to,
 *  default is 2.
 * \return returns a std::string containing converted string.
 */
std::string to_string_with_precision(const double doubleVal,
    const int precisionVal = 2);
}  // namespace util
}  // namespace glass3
#endif  // STRINGUTIL_H
