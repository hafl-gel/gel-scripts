
library(data.table)
library(ibts)

test_file <- '.test-data/20230627.TXT'

hdr <- unlist(read.table(test_file, nrows = 1, sep = ';'))
dat <- fread(cmd = paste('grep -v Date', test_file))

setnames(dat, sub('\\s+.*$', '', hdr))

# TODO: switch HIPPIE to UTC+1
dat[, c('st', 'et') := {
    Time <- fast_strptime(Date_Time, format = '%Y%m%dT%H%M%SZ', lt = FALSE)
    dT <- median(diff(Time), na.rm = TRUE)
    list(
        Time - dT,
        Time
        )
}]
dat[, c('nFlow_1', 'nFlow_2', 'Flow_1', 'Flow_2') := {
    Tair <- rowMeans(
        cbind(Temperature_outside_left, Temperature_outside_right), 
        na.rm = TRUE)
    list(
        Flow_1,
        Flow_2,
        Flow_1 * 101325 / 273.15 * (Tair + 273.15) / Air_pressure,
        Flow_2 * 101325 / 273.15 * (Tair + 273.15) / Air_pressure
    )
}]

# ibts
hippie <- as.ibts(dat, tz = 'Etc/GMT-2')

# get valve-cycle mt
valve <- dat[, .(
    mt = st[1] + (et[.N] - st[1]) / 2,
    flow1 = mean(Flow_1, na.rm = TRUE),
    flow2 = mean(Flow_2, na.rm = TRUE)
    ), by = .(Valve_Is, Cycle_Set)]

x11(width = 10)
# png('flow.png', width = 480 / 7 * 10)
plot(hippie[, 'nFlow_1'], col = 'lightgrey', ylim = c(0, 1.2e3))
lines(hippie[, 'Flow_1'])
valve[, 
    text(mt, flow1, labels = paste(Cycle_Set, Valve_Is, sep = '/'), pos = 3)
]
# dev.off()

x11(width = 10)
# png('temp.png', width = 480 / 7 * 10)
plot(hippie[, 'Temperature_inside'], col = 'orange', ylim = c(20, 55), 
    ylab = 'temperature (°C)')
lines(hippie[, 'Temperature_outside_left'], col = '#BBA529')
lines(hippie[, 'Temperature_outside_right'], col = '#BB5B29')
legend('topleft', legend = c('inside', 'outside left', 'outside right'),
    col = c('orange', '#BBA529', '#BB5B29'), lty = 1, bty = 'n')
# dev.off()


