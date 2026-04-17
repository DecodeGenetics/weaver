#pragma once
/*!
 * @file logging.hpp
 * @brief Defines the Logging class and print_xx macros.
 */

#include <cassert>
#include <functional>
#include <iomanip>
#include <memory>
#include <mutex>
#include <ostream>
#include <string_view>

#include <weaver/constants.hpp>

namespace weaver
{
/*! @brief Log message severity levels.
 *
 * @headerfile logging.hpp "weaver/logging.hpp"
 */
enum class log_severity
{
  debug,   //!< Messages are printed if `--vverbose` is specified and the code was compiled in debug mode.
  info,    //!< Messages are printed if `--verbose` or `--vverbose` is specified.
  warning, //!< Messages are always printed.
  error    //!< Messages are always printed and will abort the program.
};

//! Get a string_view containing text describing the log severity.
inline std::string_view severity2string(log_severity const severity)
{
  switch (severity)
  {
  case log_severity::debug:
    return "debug";
  case log_severity::info:
    return "info";
  case log_severity::warning:
    return "warning";
  case log_severity::error:
    return "error";
  }

  return "debug";
}

/*! @brief Stores the logging state.
 *
 * @headerfile logging.hpp "weaver/logging.hpp"
 */
class Logging
{
private:
  /*!
   * @name Types
   * @{
   */
  //! Pointer with customisable delete behaviour.
  using stream_ptr_t = std::unique_ptr<std::ostream, std::function<void(std::ostream *)>>;

  /*!
   * @}
   * @name Deleters
   * @{
   */

  //! Stream deleter that does nothing (no ownership assumed).
  static void stream_deleter_noop(std::ostream *)
  {
  }

  //! Stream deleter with default behaviour (ownership assumed).
  static void stream_deleter_default(std::ostream * ptr)
  {
    delete ptr;
  }

public:
  /*!
   * @}
   * @name Constructors
   * @{
   */

  Logging() = delete;                            //!< Deleted because of singleton pattern
  Logging(Logging const &) = delete;             //!< Deleted because of singleton pattern
  Logging(Logging &&) = delete;                  //!< Deleted because of singleton pattern
  Logging & operator=(Logging const &) = delete; //!< Deleted because of singleton pattern
  Logging & operator=(Logging &&) = delete;      //!< Deleted because of singleton pattern

  //! Construct from existing stream.
  Logging(log_severity _severity, std::ostream & _sink) noexcept;

  //! Construct from filename.
  Logging(log_severity _severity, std::string const & filename);

  /*!
   * @}
   * @name Modifying methods
   * @}
   */

  //! Call to flush the stream
  void flush_stream();

  /*!
   * @}
   * @name Public instance variables
   * @{
   */

  log_severity severity; //!< Set severity threshold for log messages.
  stream_ptr_t sink;     //!< Messages are printed into this sink.
  std::mutex mutex;      //!< Mutex that will be guarded while printing.

  /*!
   * @}
   */
};

//! Global instance of the logging state.
extern std::unique_ptr<Logging> log_singleton;

//! Prints log messages given its \a severity
template <typename... args_t>
void print_log(log_severity const severity, args_t &&... args)
{
  assert(log_singleton);
  assert(log_singleton->sink);

  if (severity < log_singleton->severity)
    return;

  // time
  auto now = std::chrono::system_clock::now();
  auto in_time_t = std::chrono::system_clock::to_time_t(now);

  // guard mutex
  {
    std::lock_guard guard{log_singleton->mutex};
    *log_singleton->sink << std::put_time(std::localtime(&in_time_t), "[%Y-%m-%d %H:%M:%S");

    auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch());
    *log_singleton->sink << '.' << std::setfill('0') << std::setw(3) << milliseconds.count() % 1000 << ']';

    // severity
    *log_singleton->sink << " <" << severity2string(severity) << "> ";

    // args
    ((*log_singleton->sink << args), ...);

    *log_singleton->sink << '\n';
#ifndef NDEBUG
    log_singleton->sink->flush();
#endif
  }
}

#ifdef NDEBUG // Release build
#  define print_debug(...) ((void)0)
#else // not NDEBUG (=> debug build)
#  define print_debug(...) print_log(weaver::log_severity::debug, __VA_ARGS__)
#endif // NDEBUG

#define print_info(...) print_log(weaver::log_severity::info, __VA_ARGS__)
#define print_warning(...) print_log(weaver::log_severity::warning, __VA_ARGS__)
#define print_error(...) print_log(weaver::log_severity::error, __VA_ARGS__)

} // namespace weaver
