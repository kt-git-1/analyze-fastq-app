import logging

from modules.logging_utils import log_tool_output, should_suppress_tool_output


def test_should_suppress_conda_version_warning():
    assert should_suppress_tool_output(
        "samtools: /home/rc/miniconda3/lib/libncursesw.so.6: "
        "no version information available (required by samtools)"
    )


def test_log_tool_output_keeps_regular_messages(caplog):
    logger = logging.getLogger("test-log-tool-output")

    with caplog.at_level(logging.DEBUG, logger="test-log-tool-output"):
        log_tool_output(logger, "samtools sort", "[bam_sort_core] merging from 4 files")

    assert "[samtools sort] [bam_sort_core] merging from 4 files" in caplog.text


def test_log_tool_output_suppresses_conda_version_warning(caplog):
    logger = logging.getLogger("test-log-tool-output-suppressed")

    with caplog.at_level(logging.DEBUG, logger="test-log-tool-output-suppressed"):
        log_tool_output(
            logger,
            "samtools sort",
            "samtools: /home/rc/miniconda3/lib/libtinfow.so.6: "
            "no version information available (required by samtools)",
        )

    assert caplog.text == ""
