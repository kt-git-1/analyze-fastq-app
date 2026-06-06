import logging
import shlex
from typing import Sequence


def format_command(cmd: Sequence[object]) -> str:
    """ログ用にコマンドを安全に1行へ整形する。"""
    return " ".join(shlex.quote(str(part)) for part in cmd)


def log_command_start(
    logger: logging.Logger,
    cmd: Sequence[object],
    step_name: str,
) -> None:
    logger.info("外部ツールを実行します: %s", step_name)
    logger.debug("実行コマンド (%s): %s", step_name, format_command(cmd))


def log_tool_output(
    logger: logging.Logger,
    step_name: str,
    line: str,
) -> None:
    logger.debug("[%s] %s", step_name, line)
