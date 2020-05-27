import os
import logging


class LoggingUtil(object):
    """
    This class encapsulates logging required by the xml to mmcif translator
    Logging can be:
        1. saved into log files with different log levels (LogFile-s)
        2. saved into streams (StreamLog-s)
        3. shown on console
    """
    class Constants(object):
        """
        Container class for the constants needed
        """
        DEFAULT_INFO_LOG_FILE_NAME = 'INFO.log'
        DEFAULT_WARN_LOG_FILE_NAME = 'WARN.log'
        DEFAULT_ERROR_LOG_FILE_NAME = 'ERROR.log'

        INFO_LOG_STRING = 'info_log_string'
        WARN_LOG_STRING = 'warn_log_string'
        ERROR_LOG_STRING = 'error_log_string'

        INFO_LOG_LEVEL = logging.INFO
        WARN_LOG_LEVEL = logging.WARNING
        ERROR_LOG_LEVEL = logging.ERROR

        def __init__(self):
            self.stream_info_const = {"name": self.INFO_LOG_STRING, "level": self.INFO_LOG_LEVEL}
            self.stream_warn_const = {"name": self.WARN_LOG_STRING, "level": self.WARN_LOG_LEVEL}
            self.stream_error_const = {"name": self.ERROR_LOG_STRING, "level": self.ERROR_LOG_LEVEL}

    class LogFile(object):
        """
        The logger can have any of three logging files with these log levels:
        1. INFO+
        2. WARN+
        3. ERROR+
        """
        def __init__(self, log_params, default_name):
            log_file = log_params.get("log_file")
            self.to_log = log_file.get("log")
            if self.to_log:
                self.filename = default_name
                name_param = log_file.get("name", None)
                if name_param:
                    self.filename = name_param
                if not os.path.exists(self.filename):
                    with open(self.filename, 'w') as f:
                        f.close()

    class StreamLog(object):
        """

        """

        def __init__(self, log_params, const):
            log_stream = log_params.get("log_stream", False)
            to_log = log_stream.get("log")
            self.hdl = None
            if to_log:
                self.level = const.get("level")
                self.open(const.get("name"))

        def open(self, stream_name):
            """

            :return:
            """
            log_string = "%(message)s"
            self.hdl = logging.StreamHandler(log_string)
            self.hdl.name = stream_name
            self.hdl.setLevel(self.level)
            formatter = logging.Formatter(log_string)
            self.hdl.setFormatter(formatter)

        def close(self):
            """

            :return:
            """
            if self.hdl:
                if not self.hdl.closed:
                    self.hdl.close()

    class Log(object):
        """
        Class for an individual log object; a list of these objects constitute logs for one entry
        """

    class EntryLogs(object):
        """
        Class that encapsulates the list of logs for an entry
        """
        _info_title = 'INFO: '
        _warn_title = 'POTENTIAL PROBLEM: '
        _error_title = 'PROBLEM: '
        _change_title = 'CHANGE MADE: '
        _not_changed_for_now_title = 'NOT CHANGED FOR NOW: '
        _validation_title = 'VALIDATION ERROR '

        def __init__(self, entry_id):
            self.id = entry_id
            self.info_logs = []
            self.warn_logs = []
            self.error_logs = []

    def __init__(self, logger_params):
        self.c = self.Constants()
        # ConsoleLogs
        self.info_log_file = None
        self.warn_log_file = None
        self.error_log_file = None
        # StreamLogs
        self.stream_info = None
        self.stream_warn = None
        self.stream_error = None
        # display logs on the console flag
        self.show_logs_on_console = False
        self.parent_logger_level = None
        self.logger = None
        if logger_params:
            self.show_on_console = logger_params.get("console")
            self.set_logging(logger_params)

    def set_logging(self, logger_params):
        """
        Sets the logger with its logging files and their permission levels

        Firstly, depending on being required, the files where the logging information is stored are created

        Secondly, the level of the parent logger is set. This value is used in
        the __del__ function in order to roll back the value. The parent logging level
        is then set to critical in order to produce the least information.


        These files can be disabled on the command line.
        :type logger_params: object
        """
        # create log file, if required
        self.info_log_file = self.LogFile(logger_params.get('info'), self.c.DEFAULT_INFO_LOG_FILE_NAME)
        self.warn_log_file = self.LogFile(logger_params.get('warn'), self.c.DEFAULT_WARN_LOG_FILE_NAME)
        self.error_log_file = self.LogFile(logger_params.get('error'), self.c.DEFAULT_ERROR_LOG_FILE_NAME)

        # note the logging level of the parent logger
        self.parent_logger_level = logging.getLogger().getEffectiveLevel()

        # keep only critical information for the translator
        logging.getLogger().setLevel(60)
        self.logger = logging.getLogger(__name__)
        if self.logger:
            # NOT SURE WHY SETTING LEVEL AGAIN
            self.logger.setLevel(logging.INFO)
            # prevent logging from being sent to the upper logger and console
            self.logger.propagate = False

            self.stream_info = self.StreamLog(logger_params.get('info'), self.c.stream_info_const)
            self.logger.addHandler(self.stream_info.hdl)
            self.stream_warn = self.StreamLog(logger_params.get('warn'), self.c.stream_warn_const)
            self.logger.addHandler(self.stream_warn.hdl)
            self.stream_error = self.StreamLog(logger_params.get('error'), self.c.stream_error_const)
            self.logger.addHandler(self.stream_error.hdl)
            return True
        else:
            return False
