
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~ function; get DOAS model info:
getDOASinfo <- function(DOASmodel, timerange = Sys.time(), tzone = "", Serial = NULL){

  timerange <- prepTimeRange(timerange,tzone)

  if (is.null(Serial)) {
    Serial <- switch(DOASmodel
      , "S1" = if(timerange[1] > parse_date_time('2021-01-01', 'Y-m-d', tz = 'Etc/GMT-1')) 'QEPB2069' else if(timerange[1] > parse_date_time("2017-01-01","Y-m-d", tz= 'Etc/GMT-1')) "QEP01232" else "AvantesS1"
      , "S2" = "QEPB0204"
      , "S3" = "QE65ProS3"
      , "S4" = "QE65ProS4"
      , "S5" = "QEPB0591"
      , "S6" = if(timerange[1] > parse_date_time("2021-01-01", "Y-m-d", tz = 'Etc/GMT-1')) 'QEPB2021' else "QEP00634"
      )    
  }

  # spectrometer information
  Spectrometer <- switch(Serial,
    # Spectrometers S1 (QEPro (neu) /Avantes (alt))
    "QEP01232" = list(
      "Spectrometer Name" = "QEPro serial-QEP01232, grating-2400"
      ,"Pixel Number" = 1044
      ,"Linearity Test Range" = c(NA,NA)
      ,"Linearity Coefficients" = c(
        "Intercept" = 0.994351,
        "Coefficient 1" = 1.12332E-7,
        "Coefficient 2" = 2.61632E-12,
        "Coefficient 3" = -1.47926E-16,
        "Coefficient 4" = 2.66477E-21,
        "Coefficient 5" = -2.31479E-26,
        "Coefficient 6" = 9.67207E-32,
        "Coefficient 7" = -1.56846E-37)
      ,"Calibration Coefficients" = c(
        136.5714762,
        0.097806116,
        -5.99541E-6,
        0)
      ),
    "AvantesS1" = list(
      "Spectrometer Name" = "AvaSpec-2048x14"
      ,"Pixel Number" = 1021
      ,"Linearity Test Range" = c(NA,NA)
      ,"Linearity Coefficients" = c(
        "Intercept" = 1,
        "Coefficient 1" = 0,
        "Coefficient 2" = 0,
        "Coefficient 3" = 0,
        "Coefficient 4" = 0,
        "Coefficient 5" = 0,
        "Coefficient 6" = 0,
        "Coefficient 7" = 0)
      ,"Calibration Coefficients" = c(
        "Intercept" = 200,
        "First Coefficient" = 0.05618,
        "Second Coefficient" = 0,
        "Third Coefficient" = 0)
      ),
    # 'neues' QE65Pro 2021
    "QEPB2069" = list(
      "Spectrometer Name" = "QE65Pro serial-QEPB2069, grating-H7"
      ,"Pixel Number" = 1044
      ,"Linearity Test Range" = c(5E3,6.5E4)
      ,"Linearity Coefficients" = c(
        "Intercept" = 0.979327,
        "Coefficient 1" = 9.93456e-7,
        "Coefficient 2" = -3.35775e-12,
        "Coefficient 3" = -1.13105e-15,
        "Coefficient 4" = 5.96016e-20,
        "Coefficient 5" = -1.56176e-24,
        "Coefficient 6" = 2.02372e-29,
        "Coefficient 7" = -1.03672e-34)
      ,"Calibration Coefficients" = c(
        "Intercept" = 137.9826798,
        "First Coefficient" = 0.096743689,
        "Second Coefficient" = -5.65527E-06,
        "Third Coefficient" = 0)
      ),
    # Spectrometer S2
    "QEPB0204" = list(
      "Spectrometer Name" = "QE65000 serial-QEPB0204, grating-2400"
      ,"Pixel Number" = 1044
      ,"Linearity Test Range" = c(NA,NA)
      ,"Linearity Coefficients" = c(
        "Intercept" = 1.03083,
        "Coefficient 1" = -3.64364E-06,
        "Coefficient 2" = 1.19736E-10,
        "Coefficient 3" = -2.0452E-15,
        "Coefficient 4" = 1.34312E-19,
        "Coefficient 5" = -5.92412E-24,
        "Coefficient 6" = 1.00883E-28,
        "Coefficient 7" = -5.87744E-34)
      ,"Calibration Coefficients" = c(
        "Intercept" = 136.7233481,
        "First Coefficient" = 0.095874504,
        "Second Coefficient" = -4.46967E-6,
        "Third Coefficient" = 0)
      ),
    # Spectrometer S3
    "QE65ProS3" = list(
      "Spectrometer Name" = "QE65000 serial-QEPB????, grating-2400"
      ,"Pixel Number" = 1044
      ,"Linearity Test Range" = c(NA,NA)
      ,"Linearity Coefficients" = c(
        "Intercept" = 0.989081,
        "Coefficient 1" = -2.03573E-06,
        "Coefficient 2" = 3.49362E-10,
        "Coefficient 3" = -2.22192E-14,
        "Coefficient 4" = 7.67257E-19,
        "Coefficient 5" = -1.49578E-23,
        "Coefficient 6" = 1.53263E-28,
        "Coefficient 7" = -6.40854E-34)
      ,"Calibration Coefficients" =c(
        "Intercept" = 124.0903261,
        "First Coefficient" = 0.146759479,
        "Second Coefficient" = -6.53796E-5,
        "Third Coefficient" = 2.31632E-8)
      ),
    # Spectrometer S4
    "QE65ProS4" = list(
      "Spectrometer Name" = "QE65000 serial-QEPB???, grating-2400"
      ,"Pixel Number" = 1044
      ,"Linearity Test Range" = c(NA,NA)
      ,"Linearity Coefficients" = c(
        "Intercept" = 0.974015,
        "Coefficient 1" = 2.70953e-06,
        "Coefficient 2" = -1.74387e-12,
        "Coefficient 3" = -1.35508e-14,
        "Coefficient 4" = 7.88328e-19,
        "Coefficient 5" = -1.94499e-23,
        "Coefficient 6" = 2.25517e-28,
        "Coefficient 7" = -1.00972e-33)
      ,"Calibration Coefficients" = c(
        "Intercept" = 120.9375051,
        "First Coefficient" = 0.158043261,
        "Second Coefficient" = -7.92905E-5,
        "Third Coefficient" = 2.86866E-8)
      ),
    # Spectrometer S5
    "QEPB0591" = list(
      "Spectrometer Name" = "QE65Pro serial-QEPB0591, grating-2400"
      ,"Pixel Number" = 1044
      ,"Linearity Test Range" = c(5E3,6.5E4)
      ,"Linearity Coefficients" = c(
        "Intercept" = 0.981007,
        "Coefficient 1" = 7.71084e-07,
        "Coefficient 2" = 9.24250e-11,
        "Coefficient 3" = -1.12364e-14,
        "Coefficient 4" = 5.36871e-19,
        "Coefficient 5" = -1.31607e-23,
        "Coefficient 6" = 1.61183e-28,
        "Coefficient 7" = -7.79958e-34)
      ,"Calibration Coefficients" =c(
        "Intercept" = 137.0514331,
        "First Coefficient" = 0.103393136,
        "Second Coefficient" = -1.64782E-5,
        "Third Coefficient" = 4.74624E-9)
      ),
    # Spectrometer S6
    # 'altes' QEPro
    "QEP00634" = list(
      "Spectrometer Name" = "QEPro serial-QEP00634, grating-????"
      ,"Pixel Number" = 1044
      ,"Linearity Test Range" = c(2E4,2E5)
      ,"Linearity Coefficients" = c(
        "Intercept" = 0.986971,
        "Coefficient 1" = 3.94404E-07,
        "Coefficient 2" = 3.09068E-12,
        "Coefficient 3" = -2.74607E-16,
        "Coefficient 4" = 4.58465E-21,
        "Coefficient 5" = -3.57549E-26,
        "Coefficient 6" = 1.35606E-31,
        "Coefficient 7" = -2.02341E-37)
      ,"Calibration Coefficients" =c(
        "Intercept" = 137.853993,
        "First Coefficient" = 0.096714797,
        "Second Coefficient" = -5.76649E-06,
        "Third Coefficient" = 0)
      ),
    # 'neues' QE65Pro 2020
    "QEPB2021" = list(
      "Spectrometer Name" = "QE65Pro serial-QEPB2021, grating-H7 Composite Blaze"
      ,"Pixel Number" = 1044
      ,"Linearity Test Range" = c(5E3,6.5E4)
      ,"Linearity Coefficients" = c(
        "Intercept" = 0.980598,
        "Coefficient 1" = 1.02522e-6,
        "Coefficient 2" = 2.01868e-11,
        "Coefficient 3" = -4.1622e-15,
        "Coefficient 4" = 1.96718e-19,
        "Coefficient 5" = -4.53588e-24,
        "Coefficient 6" = 5.17024e-29,
        "Coefficient 7" = -2.3429e-34)
      ,"Calibration Coefficients" = c(
        "Intercept" = 137.5472713,
        "First Coefficient" = 0.097457195,
        "Second Coefficient" = -6.02838E-06,
        "Third Coefficient" = 0)
      )
    )
  Spectrometer$Serial <- Serial

  # raw data structure for all models since may 2016(?)
  rawdata.structure <- list(
    "Header Lines" = 8
    ,"Time Format" = "%y%m%d%H%M%OS"
    ,"Revolver Position" = 1
    ,"Shutter Position" = 2
    ,"Acquisition start time" = 3
    ,"Acquisition stop time" = 4
    ,"TEC-Temp (degC)" = 5
    ,"Board- / Ambient-T (degC) / ambient-RH (perc)" = 6
    ,"Exposure Time (ms)" = 7
    ,"Number of accumulations" = 8
    )

  # check for old data structures
  if(DOASmodel=="S1" && timerange[1] < as.POSIXct("2016-01-01", tz=tz(timerange[1]))){
    Spectrometer$"Pixel Number" <- 2048
    rawdata.structure$"Header Lines" <- 17
    rawdata.structure$"Time Format" <- "%Y-%m-%d %H:%M:%OS"
    rawdata.structure$"Acquisition start time" <- 2
    rawdata.structure$"Acquisition stop time" <- 3  
    rawdata.structure$"Revolver Position" <- NA_integer_ 
    rawdata.structure$"Shutter Position" <- NA_integer_  
    rawdata.structure$"Exposure Time (ms)" <- 5 
  } else if(any(DOASmodel == c("S2","S3","S4"))&timerange[1] < as.POSIXct("2015-01-01", tz=tz(timerange[1]))){
    rawdata.structure$"Header Lines" <- 7
    rawdata.structure$"Board- / Ambient-T (degC) / ambient-RH (perc)" <- NA_integer_
    rawdata.structure$"Exposure Time (ms)" <- 6
    rawdata.structure$"Number of accumulations" <- 7    
  }

  # wavelength
  Spectrometer$pixel <- seq.int(Spectrometer$"Pixel Number")
  Spectrometer$wavelength <- Spectrometer$"Calibration Coefficients"[1] + Spectrometer$"Calibration Coefficients"[2]*Spectrometer$pixel + Spectrometer$"Calibration Coefficients"[3]*Spectrometer$pixel^2 + Spectrometer$"Calibration Coefficients"[4]*Spectrometer$pixel^3

  return(list(
    DOASmodel=DOASmodel
    ,Spectrometer=Spectrometer
    ,rawdata.structure=rawdata.structure
    ,timerange=timerange
    ))
}
