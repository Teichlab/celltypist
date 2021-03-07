import logging
logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)

set_level = logger.setLevel
info = logger.info
warn = logger.warning
error = logger.error
debug = logger.debug
