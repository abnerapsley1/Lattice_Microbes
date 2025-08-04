# 
# University of Illinois Open Source License
# Copyright 2008-2018 Luthey-Schulten Group,
# All rights reserved.
# 
# Developed by: Luthey-Schulten Group
#                           University of Illinois at Urbana-Champaign
#                           http://www.scs.uiuc.edu/~schulten
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the Software), to deal with 
# the Software without restriction, including without limitation the rights to 
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
# of the Software, and to permit persons to whom the Software is furnished to 
# do so, subject to the following conditions:
# 
# - Redistributions of source code must retain the above copyright notice, 
# this list of conditions and the following disclaimers.
# 
# - Redistributions in binary form must reproduce the above copyright notice, 
# this list of conditions and the following disclaimers in the documentation 
# and/or other materials provided with the distribution.
# 
# - Neither the names of the Luthey-Schulten Group, University of Illinois at
# Urbana-Champaign, nor the names of its contributors may be used to endorse or
# promote products derived from this Software without specific prior written
# permission.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL 
# THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR 
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
# OTHER DEALINGS WITH THE SOFTWARE.
# 
# Author(s): Michael J. Hallock and Joseph R. Peterson
# 
#


##############################################
### Importing Additional Necessary Modules ###
##############################################
import logging


### Create the Default LM Logger ###
# Create object LMLogger with module getLogger() method and name it 'LMLogger' #
LMLogger = logging.getLogger('LMLogger')

# Add a null handler object to LMLogger to ensure no "No handlers could be found... " warnings #
# The handlers are simply where the logging messages will be printed or saved #
# Adding a null handler means they will not be saved anywhere, but there will be no errors #
# due to the logger not being able to find a handler #
# Gracefully catch logging.NullHandler() not existing in python 2.6 and adding a new class if not available #
try:
	nullHandlerLM = logging.NullHandler()
except AttributeError:
	class NullHandler(logging.Handler):
		def emit(self, record):
			pass
	nullHandlerLM = NullHandler()
# Add the null handler to the LMLogger object #
LMLogger.addHandler(nullHandlerLM)

# Define the logging formatter (LMformatter) using logging.Formatter() method #
# Define the format to print the time, level, and message contents #
LMformatter = logging.Formatter('%(asctime)s: %(levelname)s: %(message)s')
# Set the null handler format as LMformatter #
nullHandlerLM.setFormatter(LMformatter)




### Define setLMLoggerLevel Function ###
def setLMLoggerLevel(level):
        """
        Set the level of the logger for the application

        This function enables the user to set the logger level (i.e., the category
        the user wants the message to be assigned). This is useful if you have 
	a lot of logging messages to print of the same level.

        Args:
            level (string): REQUIRED, the level the logger should report (e.g. logger.WARNING, 
	        logger.INFO, etc.)

        Returns:
            none

        """
	
        LMLogger.setLevel(level)





def setLMLogFile(filename, level=logging.DEBUG):
        """
        Set the file name for logging the information

        This function enables the user to set the logger file which will
	specify the file name where all logging information is printed/saved.

        Args:
            filename (string): REQUIRED, the name of the file for the logging information
	        to be stored.
	    level (string): The level for the logger to set the messages at. This needs
                to be in the following format "logging.LEVEL", where LEVEL can be DEBUG,
		INFO, WARNING, ERROR, CRITICAL, etc.

        Returns:
            none

        """

	
	fileH = logging.FileHandler(filename, mode='w')
	fileH.setLevel(level)
	fileH.setFormatter(LMformatter)
	LMLogger.removeHandler(nullHandlerLM)
	LMLogger.addHandler(fileH)





def setLMLogConsole(level=logging.DEBUG):
	"""Set the logger to write to the console as the code is working

    Args:
        level:
            The level of information to log
    """
	consoleH = logging.StreamHandler() # Defaults to sys.stderr
	consoleH.setLevel(level)
	consoleH.setFormatter(LMformatter)
	LMLogger.removeHandler(nullHandlerLM)
	LMLogger.addHandler(consoleH)

