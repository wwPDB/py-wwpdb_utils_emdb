
import pytest
import sys
from translator_classes.EMDBXmlToCifTranslator import EMDBXmlToCifTranslator
from translator_classes.CIF import CIF
from translator_classes.EMDBMetadata import EMDBMetadata
from translator_classes.LoggingUtil import LoggingUtil

builtins_str = None
file_nf_error = Exception
if sys.version_info.major == 3:
    from unittest.mock import MagicMock
    builtins_str = "builtins"
    file_nf_error = FileNotFoundError
elif sys.version_info.major == 2:
    from mock import MagicMock
    builtins_str = "__builtin__"
    file_nf_error = IOError

test_logger_params = {
    "info": {
        "log_file": {
            "log": True,
            "name": "/Users/sanja/IdeaProjects/emdb-xml2cif-translator/logs/info.log"
        },
        "log_stream": {
            "log": True,
            "name": None
        }
    },
    "warn": {
        "log_file": {
            "log": True,
            "name": "/Users/sanja/IdeaProjects/emdb-xml2cif-translator/logs/warn.log"
        },
        "log_stream": {
            "log": True,
            "name": None
        }
    },
    "error": {
        "log_file": {
            "log": True,
            "name": "/Users/sanja/IdeaProjects/emdb-xml2cif-translator/logs/error.log"
        },
        "log_stream": {
            "log": True,
            "name": None
        }
    },
    "console": True
}

input_file = 'emd-00000-v30.xml'
output_file = 'emd-00000.cif'
emdb_id = '00000'


@pytest.fixture()
def translator():
    the_translator = EMDBXmlToCifTranslator(test_logger_params)
    return the_translator


@pytest.fixture()
def emdb_metadata():
    the_metadata = EMDBMetadata()
    the_metadata.cif = CIF(output_file)
    return the_metadata


@pytest.fixture()
def mock_open(monkeypatch):
    mock_file = MagicMock()
    mock_open = MagicMock(return_value=mock_file)
    open_str = builtins_str + ".open"
    monkeypatch.setattr(open_str, mock_open)
    return mock_open


def test_can_read_input_file(monkeypatch, translator, emdb_metadata):
    mock_exists = MagicMock(return_value=True)
    monkeypatch.setattr("os.path.exists", mock_exists)
    mock_xml_tree = MagicMock()
    emdb_metadata.xml_tree_created = False
    translator.emdb_data = emdb_metadata
    monkeypatch.setattr(emdb_metadata, 'parse_into_xml_tree', mock_exists)
    translator.emdb_data.xml_tree = mock_xml_tree
    result = translator.read_emdb_header_file(input_file)
    assert result is True, 'Input file' + input_file + ' cannot be read'


def test_read_input_file_exception(monkeypatch, translator, emdb_metadata):
    mock_exists = MagicMock(return_value=False)
    monkeypatch.setattr("os.path.exists", mock_exists)
    with pytest.raises(Exception) as e:
        result = translator.read_emdb_header_file(input_file, emdb_metadata)
        assert e == file_nf_error, 'Reading the input XML file raises an exception but it is not ' + str(e)
        assert result is False, 'An exception is raised but the return is set to success'


def test_can_write_cif_file(translator, emdb_metadata):
    translator._EMDBXmlToCifTranslator__emdb_header_file_translated = True
    translator.emdb_data = emdb_metadata
    translator.emdb_data.cif.filename = input_file
    result = translator.write_mmcif_file()
    assert result is True, 'Output mmcif file is not written'


def test_translate_all_in(monkeypatch, translator):
    mock_exists = MagicMock(return_value=True)
    monkeypatch.setattr("os.path.exists", mock_exists)
    mock_return = MagicMock(return_value=True)
    translator.read_emdb_header_file = mock_return
    translator.write_mmcif_file = mock_return
    translator.translate_xml_to_cif(input_file, output_file)


def test_translate_with_only_input_file(translator):
    mock_return = MagicMock(return_value=True)
    translator.read_emdb_header_file = mock_return
    translator.translate_xml_to_cif(input_file)


def test_can_set_up_a_logger(translator):
    translator.logging = LoggingUtil(test_logger_params)
    assert translator.logging is not None
    assert translator.logging.show_on_console is test_logger_params.get('console')

    info_params = test_logger_params.get('info')
    warn_params = test_logger_params.get('warn')
    error_params = test_logger_params.get('error')

    if info_params:
        assert translator.logging.info_log_file.to_log is info_params.get('log_file').get("log")
        assert translator.logging.info_log_file.filename is not None
    if warn_params:
        assert translator.logging.warn_log_file.to_log is warn_params.get('log_file').get("log")
        assert translator.logging.warn_log_file is not None
    if error_params:
        assert translator.logging.error_log_file.to_log is error_params.get('log_file').get("log")
        assert translator.logging.error_log_file is not None


def test_can_set_up_translation_log(translator):
    translator.translation_log = LoggingUtil.EntryLogs(emdb_id)
    assert translator.translation_log.id is emdb_id
    assert translator.translation_log.info_logs == []
    assert translator.translation_log.warn_logs == []
    assert translator.translation_log.error_logs == []


# def test_can_output_translation_logs(translator):
#     translator.logging = LoggingUtil(test_logger_params)
#     translator.logging.translation_log = LoggingUtil.EntryLogs(emdb_id)
#     assert translator.output_translation_log() is True



