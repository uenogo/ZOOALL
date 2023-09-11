import sys, os
import socket

sys.path.append("/isilon/BL32XU/BLsoft/PPPP/10.Zoo/Libs")

import BeamsizeConfig, Date
import logging
import logging.config
import Zoo

if __name__ == "__main__":
    config_dir = "/isilon/blconfig/bl32xu/"

    # kuntaro_log
    d = Date.Date()
    time_str = d.getNowMyFormat(option="date")
    logname = "/isilon/BL32XU/BLsoft/PPPP/10.Zoo/ZooLogs/debug_%s.log" % time_str
    logging.config.fileConfig('/isilon/BL32XU/BLsoft/PPPP/10.Zoo/Libs/logging.conf', defaults={'logfile_name': logname})
    logger = logging.getLogger('ZOO')
 
    zoo = Zoo.Zoo()
    zoo.connect()

    while(True):
        wavelength = zoo.getWavelength()
        print(("wavelength:", wavelength))
        beamsize_index = zoo.getBeamsize()
        print(("beamsize_index:", beamsize_index))
        zoo.waitSPACE()
