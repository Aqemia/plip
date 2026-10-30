import inspect
import logging

# Quiet by default when used as a library; plipcmd overrides this from the CLI flags.
logging.getLogger('plip').setLevel(logging.WARNING)


def get_logger():
    """
    Configures a base logger and returns a module-specific sub-logger of the calling module.
    """
    frame = inspect.stack()[1]
    module_name = inspect.getmodule(frame[0]).__name__
    if module_name != '__main__':
        logger = logging.getLogger(module_name)
        root = logging.getLogger()
        if not root.handlers:
            ch = logging.StreamHandler()
            formatter = logging.Formatter('%(asctime)s [%(levelname)s] [%(filename)s:%(lineno)d] %(name)s: %(message)s')
            ch.setFormatter(formatter)
            root.addHandler(ch)
    else:
        logger = logging.getLogger('plip')

    return logger
