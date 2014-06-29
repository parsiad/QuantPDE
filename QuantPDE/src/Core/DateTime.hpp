#ifndef QUANT_PDE_CORE_DATETIME
#define QUANT_PDE_CORE_DATETIME

#include <cstring>  // std::strlen
#include <ctime>    // std::time_t, std::tm, gmtime_r, timegm, asctime_r
#include <string>   // std::string
#include <iostream> // std::ostream

namespace QuantPDE {

/**
 * A convenience class for manipulating dates.
 */
class DateTime final {

	std::time_t time;
	std::tm details;

public:

	enum Weekday {
		Sunday    = 1,
		Monday    = 2,
		Tuesday   = 3,
		Wednesday = 4,
		Thursday  = 5,
		Friday    = 6,
		Saturday  = 7
	};

	enum Month {
		January   = 1,
		February  = 2,
		March     = 3,
		April     = 4,
		May       = 5,
		June      = 6,
		July      = 7,
		August    = 8,
		September = 9,
		October   = 10,
		November  = 11,
		December  = 12
	};

	typedef int Seconds;

	typedef int Minutes;

	typedef int Hours;

	typedef int Day;

	typedef int Year;

	typedef int YearDay;

	/**
	 * Initialize using a UNIX timestamp.
	 * @param time A UNIX timestamp.
	 */
	DateTime(std::time_t time = 0) noexcept : time(time) {
		gmtime_r(&time, &details);
	}

	/**
	 * Initialize using a string date.
	 * @param str The date-time string.
	 * @param format The format used to perform the parsing. See <ctime>.
	 */
	DateTime(const std::string &str,
			const std::string &format = "%Y-%m-%d %T") noexcept {
		strptime(str.c_str(), format.c_str(), &details);
		time = timegm(&details); //mktime(&details);
	}

	/**
	 * Constructor.
	 */
	DateTime(int year, int month, int day, int hours = 0, int minutes = 0,
			int seconds = 0) noexcept {
		details.tm_sec = seconds;
		details.tm_min = minutes;
		details.tm_hour = hours;
		details.tm_mday = day;
		details.tm_mon = month - 1;
		details.tm_year = year - 1900;

		time = timegm(&details);
		gmtime_r(&time, &details);
	}

	/**
	 * Copy constructor.
	 */
	DateTime(const DateTime &that) noexcept : time(that.time),
			details(that.details) {
	}

	/**
	 * Copy assignment operator.
	 */
	DateTime &operator=(const DateTime &that) & noexcept {
		time = that.time;
		details = that.details;
		return *this;
	}

	/**
	 * @return The UNIX timestamp.
	 */
	std::time_t timestamp() const {
		return time;
	}

	/**
	 * @return Seconds (0-60).
	 */
	Seconds seconds() const {
		return details.tm_sec;
	}

	/**
	 * @return Minutes (0-59).
	 */
	Minutes minutes() const {
		return details.tm_min;
	}

	/**
	 * @return Hours (0-23).
	 */
	Hours hours() const {
		return details.tm_hour;
	}

	/**
	 * @return Day of the month (1-31).
	 */
	Day day() const {
		return details.tm_mday;
	}

	/**
	 * @return Month.
	 */
	Month month() const {
		return Month(details.tm_mon + 1);
	}

	/**
	 * @return Year.
	 */
	Year year() const {
		return 1900 + details.tm_year;
	}

	/**
	 * @return Day of the week.
	 */
	Weekday weekday() const {
		return Weekday(details.tm_wday + 1);
	}

	/**
	 * @return Day in the year (0-365, 1 Jan = 0).
	 */
	YearDay yearDay() const {
		return details.tm_yday;
	}

	/**
	 * @return True if and only if this date is in daylight savings.
	 */
	bool daylightSavings() const {
		return details.tm_isdst;
	}

	/**
	 * Prettifies and prints the date to an output stream.
	 * @param os The output stream.
	 * @param axis The axis.
	 */
	friend std::ostream &operator<<(std::ostream &os,
			const DateTime &dateTime) {
		// Create string
		char buffer[26];
		asctime_r(&dateTime.details, buffer);

		// Remove new line
		buffer[std::strlen(buffer) - 1] = '\0';

		// Print and return
		os << buffer;
		return os;
	}

};

/**
 * @return The difference between the two times, in years.
 */
Real operator-(const DateTime &a, const DateTime &b) {
	return (double) (a.timestamp() - b.timestamp()) / 60. / 60. / 24.
			/ 365. ;
}

/**
 * @return True if and only if the two objects point to the same date-time.
 */
bool operator==(const DateTime &a, const DateTime &b) {
	return a.timestamp() == b.timestamp();
}

/**
 * @see QuantPDE::DateTime::operator==
 */
bool operator!=(const DateTime &a, const DateTime &b) {
	return a.timestamp() != b.timestamp();
}

/**
 * @return True if and only if the first argument points to an earlier date than
 *         the second.
 */
bool operator<(const DateTime &a, const DateTime &b) {
	return a.timestamp() < b.timestamp();
}

/**
 * @see QuantPDE::DateTime::operator<
 * @see QuantPDE::DateTime::operator==
 */
bool operator>(const DateTime &a, const DateTime &b) {
	return a.timestamp() > b.timestamp();
}

/**
 * @see QuantPDE::DateTime::operator<
 * @see QuantPDE::DateTime::operator==
 */
bool operator<=(const DateTime &a, const DateTime &b) {
	return a.timestamp() <= b.timestamp();
}

/**
 * @see QuantPDE::DateTime::operator<
 * @see QuantPDE::DateTime::operator==
 */
bool operator>=(const DateTime &a, const DateTime &b) {
	return a.timestamp() >= b.timestamp();
}

}

#endif

