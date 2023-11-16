
# NOTES:

# TODO:

## 0. header ----------------------------------------

## 1. functions ----------------------------------------

## 2. testing ----------------------------------------

path_lfe <- '~/LFE/02_Daten/9-Messventilator'
fan_id <- '1'
path_prep <- file.path(path_lfe, 'messventilator-%s0mm-%s/frequi-%s_%s_%i_%i_%i.csv')

# get diameter
fan_dia <- switch(as.character(fan_id)
    , '1' = 63
    , '2' = 
    , '3' = 
    , '4' = 
    , '5' = 82
    , 92
)

# create path
path_fan
'/messventilator-630mm-1/frequi-1_63_2023_11_14.csv'
