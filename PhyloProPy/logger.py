import logging

def phyloprofile_logger(log_file='./phyloprofile.log', debug=False, silent=False):
    
    # parse logging level
    if debug:
        level = logging.DEBUG
    elif silent:
        level = logging.WARNING
    else:
        level = logging.INFO
    
    # Create a logger instance
    logger = logging.getLogger('phyloprofile')
    if logger.hasHandlers():
        logger.handlers.clear()
    logger.setLevel(level)
    
    # Create a file handler and a console handler
    #file_handler = logging.FileHandler(log_file)
    console_handler = logging.StreamHandler()
    
    # Create a formatter
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    
    # Set the formatter for both handlers
    #file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)
    
    # Add both handlers to the logger
    #logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    return logger