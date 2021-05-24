library("DescTools")
library("fields")
library("ggplot2")
library("lubridate")
library("raster")
library("RColorBrewer")
library("zoo")

tab = read.csv("Data_LUHB-ULB_200521.csv", head=T, sep=",")

# 1. Estimating the median and mean Ct values through time

bootstrapsCI = function(vS, type="median", conf.level=0.90, R=999)
	{
		differences = rep(NA, R); vS = vS[which(!is.na(vS))]
		for (i in 1:R)
			{
				if (type == "median") differences[i] = median(sample(vS,length(vS),replace=T))-median(vS)
				if (type == "mean") differences[i] = mean(sample(vS,length(vS),replace=T))-median(vS)
			}
		return(median(vS)+quantile(differences, probs=c(0.5-(conf.level/2),0.5+(conf.level/2))))
	}

confidence_interval = FALSE; confidence_interval = TRUE
dates = dmy(tab[,"Date.prelv"]); years = c(rep(2020,10), rep(2021,5))
months = c("03","04","05","06","07","08","09","10","11","12","01","02","03","04","05")
startDays = c(1, rep(1,length(months)-1))
endDays = c(31,30,31,30,31,31,30,31,30,31,31,28,31,30,15)
labels1 = c(); labels2 = c(); labels3 = c()
cases1 = list(); cases2 = list(); days = list(); n = 0
cts1_a = list(); cts1_a_low = list(); cts1_a_high = list()
cts1_b = list(); cts1_b_low = list(); cts1_b_high = list()
cts2_a = list(); cts2_a_low = list(); cts2_a_high = list()
cts2_b = list(); cts2_b_low = list(); cts2_b_high = list()
for (i in 1:length(months))
	{
		for (j in startDays[i]:endDays[i])
			{
				n = n+1
				if (j == 1)
					{
						labels1 = c(labels1, n)
						labels2 = c(labels2, paste("01",months[i],years[i],sep="/"))
					}
				if (j < 10) j = paste0("0",j)
				labels3 = c(labels3, paste(j,months[i],years[i],sep="/"))
				date = dmy(paste(j,months[i],years[i],sep="/"))
				sub = tab[which(dates==date),]; days[[n]] = date
				indices = which(sub[,"RSLT"]=="POS"); cases1[[n]] = length(indices)
				indices = which((sub[,"RSLT"]=="POS")&(grepl("M 2000 Real Time",sub[,"Analyseur"])))
					# the other analysors did not lead to comparable Ct values
				vS = sub[indices,"Ct"]
				cts1_a[[n]] = median(vS, na.rm=T)
				if (length(vS) > 0)
					{
						if (confidence_interval == TRUE)
							{
								if (sum(!is.na(vS)) > 1)
									{
										ci = MedianCI(vS, conf.level=0.90, sides="two.sided", na.rm=T, method="boot", R=999)
										cts1_a_low[[n]] = ci["lwr.ci"]; cts1_a_high[[n]] = ci["upr.ci"]
										ci = bootstrapsCI(vS, type="median", conf.level=0.90, R=999)
										cts1_a_low[[n]] = ci[1]; cts1_a_high[[n]] = ci[2]
									}	else		{
										cts1_a_low[[n]] = NA; cts1_a_high[[n]] = NA
									}
							}	else	{		
								if (sum(!is.na(cts1_a[[n]])) > 0)
									{
										cts1_a_low[[n]] = quantile(vS, 0.25, na.rm=T)
										cts1_a_high[[n]] = quantile(vS, 0.75, na.rm=T)
									}	else	{
										cts1_a_low[[n]] = NA; cts1_a_high[[n]] = NA
									}
							}
					}	else		{
						cts1_a_low[[n]] = NA; cts1_a_high[[n]] = NA
					}
				cts1_b[[n]] = mean(vS, na.rm=T)
				if (length(vS) > 0)
					{
						if (confidence_interval == TRUE)
							{
								if (sum(!is.na(vS)) > 1)
									{
										ci = MeanCI(vS, conf.level=0.90, sides="two.sided", na.rm=T, method="boot", R=999)
										cts1_b_low[[n]] = ci["lwr.ci"]; cts1_b_high[[n]] = ci["upr.ci"]
										ci = bootstrapsCI(vS, type="mean", conf.level=0.90, R=999)
										cts1_b_low[[n]] = ci[1]; cts1_b_high[[n]] = ci[2]
									}	else		{
										cts1_b_low[[n]] = NA; cts1_b_high[[n]] = NA
									}
							}	else	{
								if (sum(!is.na(cts1_a[[n]])) > 0)
									{
										cts1_b_low[[n]] = quantile(vS, 0.25, na.rm=T)
										cts1_b_high[[n]] = quantile(vS, 0.75, na.rm=T)
									}	else	{
										cts1_b_low[[n]] = NA; cts1_b_high[[n]] = NA
									}
							}
					}	else		{
						cts1_b_low[[n]] = NA; cts1_b_high[[n]] = NA
					}
			}
	}

# 2. Co-plotting the mean Ct values on GEES phase diagrams

slidingWindow = 14; gap_Ct = 0
for (i in (slidingWindow+1):length(cases1))
	{
		cases2[[i]] = mean(unlist(cases1[(i-slidingWindow):i]), na.rm=T)
		cts2_a[[i]] = -mean(unlist(cts1_a[(i-slidingWindow):i]), na.rm=T)
		cts2_a_low[[i]] = -mean(unlist(cts1_a_low[(i-slidingWindow):i]), na.rm=T)
		cts2_a_high[[i]] = -mean(unlist(cts1_a_high[(i-slidingWindow):i]), na.rm=T)
		cts2_b[[i]] = -mean(unlist(cts1_b[(i-slidingWindow):i]), na.rm=T)
		cts2_b_low[[i]] = -mean(unlist(cts1_b_low[(i-slidingWindow):i]), na.rm=T)
		cts2_b_high[[i]] = -mean(unlist(cts1_b_high[(i-slidingWindow):i]), na.rm=T)
	}
for (i in 1:length(cts2_a))
	{
		if (is.null(cts2_a[[i]])) cts2_a[[i]] = NA
		if (is.null(cts2_b[[i]])) cts2_b[[i]] = NA
	}

date_ranges = list() # 1° division by epidemic phases:
date_ranges[[1]] = c("2020-03-01","2020-06-30")
date_ranges[[2]] = c("2020-07-01","2020-08-31")
date_ranges[[3]] = c("2020-09-01","2021-01-31")
date_ranges[[4]] = c("2021-02-01","2021-04-30")
date_ranges = list() # 2° division by epidemic phases:
date_ranges[[1]] = c("2020-03-01","2020-07-31")
date_ranges[[2]] = c("2020-07-14","2020-09-01")
date_ranges[[3]] = c("2020-09-01","2021-02-01")
date_ranges[[4]] = c("2021-02-01","2021-05-15")
load("Background_C_Faes.RData")
weekend_correction = TRUE
plottingDashedLines = FALSE
usingMedianCtValues = FALSE
usingMedianCtValues = TRUE

if (usingMedianCtValues == TRUE) vS = -unlist(cts2_a)
if (usingMedianCtValues != TRUE) vS = -unlist(cts2_b)
last_colour = colorRampPalette(brewer.pal(9,"Blues"))(66)[66]
colour_scale = rev(c(colorRampPalette(brewer.pal(9,"Blues"))(66)[1:66],rep(last_colour,34)))
colour_scale = rev(colorRampPalette(brewer.pal(9,"Blues"))(101)[1:101]) # original colour scale
cols1 = colour_scale[(((vS-min(vS,na.rm=T))/(max(vS,na.rm=T)-min(vS,na.rm=T)))*100)+1]
if (gap_Ct > 0) cols1 = c(rep(NA,gap_Ct), cols1) # to manually add a gap of minus 14 days on the Ct values
if (usingMedianCtValues == TRUE) pdf("Figure_3_NEW.pdf", width=10, height=8)
if (usingMedianCtValues != TRUE) pdf("Figure_S2_NEW.pdf", width=10, height=8)
par(mfrow=c(2,2), oma=c(1,1,0,0), mar=c(2.2,3.1,1,1), lwd=0.2, col="gray30")
for (h in 1:length(date_ranges))
	{
		hosp = 10^seq(0.9,3,0.0025); gr = seq(0.85,1.25,0.0025)
		breaks = c(0,0.025*2000,0.15*2000,0.25*2000,0.50*2000,2000,10^10)
		breaks = c(0,50,303,528,987,1502,10^10)
		colours = c(hcl.colors(6,"Heat",rev=T)[1:5], "#aa3c3c")
		image(log10(hosp), gr, matrix(result[,3],nrow=length(hosp),byrow=T), axes=F, ann=F, col=colours, breaks=breaks)
		xaxis_seq = c(0,10,20,40,75,150,300,600); yaxis_seq = c(0.9,1,1.025,1.050,1.1,1.2)
		axis(1, at=log10(xaxis_seq), labels=xaxis_seq, cex.axis=0.75, lwd=0, lwd.tick=0.2, tck=-0.015, mgp=c(0,0.15,0), col="gray30", col.tick="gray30", col.axis="gray30")
		axis(2, at=yaxis_seq, labels=c("-10%","0.0%","2.5%","5.0%","10%","20%"), las=1, cex.axis=0.75, lwd=0, lwd.tick=0.2, 
			 tck=-0.015, mgp=c(0,0.5,0), col="gray30", col.tick="gray30", col.axis="gray30")
		title(xlab="New hospitalisations", cex.lab=0.9, mgp=c(1.1,0,0), col.lab="gray30")
		title(ylab="New hospitalisations daily ratio", cex.lab=0.9, mgp=c(2.3,0,0), col.lab="gray30")
		mtext(date_ranges[[h]][1], line=-1.2, at=1.09, cex=0.7, col="gray30")
		mtext(date_ranges[[h]][2], line=-2, at=1.09, cex=0.7, col="gray30")
		if (plottingDashedLines == TRUE)
			{
				segments(0, 1, 1000, 1, lty="dashed", lwd=2)
				segments(0, 1.025, 1000, 0.95, lty="dotted", lwd=1)
				segments(0, 1.05, 1000, 0.95, lty="dotted", lwd=1)
				segments(log10(75), 0.8, log10(75), 1.5, lty="dotted", lwd=2)
				segments(log10(150), 0.8, log10(150), 1.5, lty="dotted", lwd=2)
			}
		box(lwd=0.2, col="gray30")
		contour(log10(hosp), gr, matrix(result[,3],nrow=length(hosp),byrow=T), add=T, cex=1.0, col="gray10",
				levels=breaks, labels=c("0","50 ICU beds","Phase 1A","Phase 1B","Phase 2A","Phase 2B"), lwd=0.3)
		dta = read.csv("https://epistat.sciensano.be/Data/COVID19BE_hosp.csv")
		dta_agg1 = aggregate(NEW_IN ~ DATE, dta, sum)
		dta_agg2 = aggregate(TOTAL_IN_ICU ~ DATE, dta, sum)
		dta_agg = merge(dta_agg1, dta_agg2, by="DATE")
		dta_agg$DATE = as.character(dta_agg$DATE)
		N = nrow(dta_agg); dta_agg$ICU_IN2WEEKS = NA
		dta_agg$ICU_IN2WEEKS[1:(N-14)] = dta_agg$TOTAL_IN_ICU[15:N]
		dta_agg$new_in_mean = rollmean(dta_agg$NEW_IN, 7, align="right", fill=NA)
		windows = 14 # using 14 instead of 12 to remove weekend effects
		dta_agg$new_in_growth = rep(NA,nrow(dta_agg))
		for (i in windows:nrow(dta_agg))
			{
				temp_df = as.data.frame(cbind(0:(windows-1),dta_agg$NEW_IN[(i-windows+1):i]))
				names(temp_df) = c("day","new_in")
				dates = dta_agg$DATE[(i-windows+1):i]
				if (weekend_correction == TRUE)
					{
						temp_df$weekend = as.numeric(weekdays(as.Date(dates))%in%c("Monday","Sunday"))
						holidays = c("2020-04-13","2020-05-01","2020-05-21","2020-07-21","2020-11-11","2020-12-25","2021-01-01","2021-04-05")
						temp_df$weekend[dates%in%holidays] = 1
						mylm = lm(log10(new_in) ~ day + weekend, temp_df)
					}	else	{
						mylm = lm(log10(new_in) ~ day, temp_df)
					}
				dta_agg$new_in_growth[i] = 10^coefficients(mylm)[[2]]
			}
		if (!is.na(date_ranges[[h]][1])) dta_agg = subset(dta_agg,DATE>=as.Date(date_ranges[[h]][1]))
		if (!is.na(date_ranges[[h]][2])) dta_agg = subset(dta_agg,DATE<=as.Date(date_ranges[[h]][2]))
		dta_agg_s = dta_agg[!is.na(dta_agg$new_in_growth),]
		dta_agg_s$nday = 1:(nrow(dta_agg_s))
		dta_agg_s$transparency = dta_agg_s$nday/nrow(dta_agg_s)/2+0.5
		dta_agg_s$size = dta_agg_s$TOTAL_IN_ICU/500
		xaxis_seq = c(12.5,25,50,100,200,400,600)
		yaxis_seq = seq(0.9,1.30,0.05); pos.date = c(1,1); cols2 = rep(NA, dim(dta_agg)[1])
		for (i in 1:dim(dta_agg_s)[1])
			{
				date = unlist(strsplit(as.character(dta_agg_s[i,"DATE"]),"-"))
				date = paste(date[3],date[2],date[1],sep="/")
				if (date%in%labels3) cols2[i] = cols1[which(labels3==date)]
			}
		lines(log10(dta_agg_s$new_in_mean), dta_agg_s$new_in_growth, lty=1, col="gray50", lwd=1)
		for (i in 1:length(dta_agg_s$new_in_mean))
			{
				if (!is.na(cols2[i]))
					{
						points(log10(dta_agg_s$new_in_mean[i]), dta_agg_s$new_in_growth[i], pch=16, col=cols2[i], cex=1.0)
						points(log10(dta_agg_s$new_in_mean[i]), dta_agg_s$new_in_growth[i], pch=1, col="gray30", cex=1.0)		
					}	else	{
						points(log10(dta_agg_s$new_in_mean[i]), dta_agg_s$new_in_growth[i], pch=16, col="gray50", cex=0.6)
					}
			}
		text(log10(dta_agg_s$new_in_mean[1]), dta_agg_s$new_in_growth[1], format(as.Date(dta_agg_s$DATE[1]),"%d/%m"), pos=pos.date[1], cex=0.8, col="gray30")
		text(log10(dta_agg_s$new_in_mean[nrow(dta_agg_s)]), dta_agg_s$new_in_growth[nrow(dta_agg_s)], format(as.Date(dta_agg_s$DATE[nrow(dta_agg_s)]),"%d/%m"), 
			 pos=pos.date[2], cex=0.8, col="gray30")
		if (h == 2)
			{
				legendRast = raster(as.matrix(c(-max(vS,na.rm=T),-min(vS,na.rm=T))))
				plot(legendRast, legend.only=T, col=rev(colour_scale), legend.width=0.5, legend.shrink=0.3,
					 smallplot=c(0.88,0.89,0.60,0.92), alpha=1, horizontal=F, legend.args=list(text="", cex=0.7, line=0.5, col="gray10"),
					 axis.args=list(cex.axis=0.8, lwd=0, lwd.tick=0.3, tck=-1.2, col.axis="gray10", col.tick="gray10", line=0, mgp=c(0,0.5,0),
					 at=seq(-14,-24,-2), labels=seq(14,24,2)))
			}
	}
dev.off()

# 4. Plotting the median and mean Ct values through time

startingDay = 61; xS = startingDay:length(cases1); cols1 = list(); cols2 = list()
cols1[[1]] = rgb(150,150,150,255,maxColorValue=255); cols2[[1]] = rgb(150,150,150,120,maxColorValue=255)
cols1[[2]] = rgb(222,67,39,255,maxColorValue=255); cols2[[2]] = rgb(222,67,39,150,maxColorValue=255); cols1[[3]] = "green3"
events_1 = decimal_date(dmy(c("10-05-2020","11-05-2020","18-05-2020","08-06-2020","01-07-2020","03-07-2020","29-07-2020",
							  "01-09-2020","09-10-2020","19-10-2020","02-11-2020","16-11-2020","01-12-2020","04-01-2021",
							  "27-01-2021","01-03-2021","08-03-2021","27-03-2021","19-04-2021","26-04-2021")))
	# 10/05 - gathering of maximum 4 people
	# 11/05 - re-opening of non-essential shops
	# 18-05 - progressive re-opening of schools
	# 08-06 - one-day trips in Belgium re-allowed
	# 08-06 - re-opening of bars and restaurants,
	# 		  gathering of max. 10 persons per week
	# 01-07 - re-opening of theme parks, spas,
	#		  restaurants, and bars
	# 03-07 - gathering of max. 15 persons per week
	# 29-07 - max. 5 close contacts/family/week,
	# 		  home working recommended
	# 01-09 - schools re-opening
	# 09-10 - bars and restaurant closed at 11 pm
	# 19-10 - closing of bars/restaurants + curfew
	# 02-11 - partial lockdown (1 close contact/pers.)
	# 16-11 - schools re-opening (after two weeks)
	# 01-12 - re-opening of non-essential shops
	# 04-01 - schools re-opening (after two weeks)
	# 27-01 - non-essential trips abroad forbidden 
	# 01-03 - re-opening of contact professions
	# 08-03 - gathering of max. 10 persons outside
	# 27-03 - closing of schools, non-essential shops,
	# 		  contact professions, max. 4 pers. outside
	# 19-04 - schools re-opening (after 3 weeks)
	# 26-04 - re-opening of contact professions,
	# 		  non-essential shops by appointment
events_2 = rep(NA, length(events_1))
for (i in 1:length(events_1))
	{
		events_2[i] = which(decimal_date(dmy(gsub("\\/","-",labels3)))==events_1[i])
	}

pdf("Figure_2_NEW.pdf", width=9, height=6); par(mfrow=c(2,1), oma=c(0,0,0,0), mar=c(2.0,1.8,0.2,1.5), col="gray30", lwd=0.2)
slidingWindow = 14
for (i in (slidingWindow+1):length(cases1))
	{
		cases2[[i]] = mean(unlist(cases1[(i-slidingWindow):i]), na.rm=T)
		cts2_a[[i]] = -mean(unlist(cts1_a[(i-slidingWindow):i]), na.rm=T)
		cts2_a_low[[i]] = -mean(unlist(cts1_a_low[(i-slidingWindow):i]), na.rm=T)
		cts2_a_high[[i]] = -mean(unlist(cts1_a_high[(i-slidingWindow):i]), na.rm=T)
		cts2_b[[i]] = -mean(unlist(cts1_b[(i-slidingWindow):i]), na.rm=T)
		cts2_b_low[[i]] = -mean(unlist(cts1_b_low[(i-slidingWindow):i]), na.rm=T)
		cts2_b_high[[i]] = -mean(unlist(cts1_b_high[(i-slidingWindow):i]), na.rm=T)
	}
plot(xS, unlist(cases2[startingDay:length(cases1)]), col=NA, lwd=0.5, type="l", axes=F, ann=F, xlim=c(min(xS),max(xS)))
yy_l = c(rep(0,length(xS)),rev(unlist(cases2[startingDay:length(cases1)])))
xx_l = c(xS,rev(xS)); polygon(xx_l, yy_l, col=cols2[[1]], border=0)
for (i in 1:length(events_2)) abline(v=events_2[i], col="gray50", lwd=0.7, lty=2)
axis(side=1, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.02, col.axis="gray30", mgp=c(0,0.03,0), at=labels1, labels=labels2)
axis(side=2, lwd.tick=0.2, cex.axis=0.7, lwd=0.2, tck=-0.02, col.axis="gray30", mgp=c(0,0.20,0), at=seq(0,200,50), pos=labels1[3]-3)
par(new=T); yMin = min(c(unlist(cts2_a_low),unlist(cts2_a_high)), na.rm=T); yMax = max(c(unlist(cts2_a_low),unlist(cts2_a_high)), na.rm=T)
plot(xS, unlist(cts2_a[startingDay:length(cases1)]), col=cols1[[2]], lwd=1.5, type="l", axes=F, ann=F, xlim=c(min(xS),max(xS)), ylim=c(-30,-8))
yy_l = c(rep(unlist(cts2_a_high[startingDay:length(cases1)])),rev(unlist(cts2_a_low[startingDay:length(cases1)])))
xx_l = c(xS,rev(xS)); polygon(xx_l, yy_l, col=cols2[[2]], border=0)
axis(side=4, lwd.tick=0.2, cex.axis=0.7, lwd=0.2, tck=-0.02, col.axis=cols1[[1]], mgp=c(0,0.10,0), at=seq(-40,0,5), labels=seq(40,0,-5))
lines(xS, unlist(cts2_b[startingDay:length(cases1)]), col=cols1[[2]], lwd=1.5, lty=2)
lines(xS, unlist(cts2_b_low[startingDay:length(cases1)]), col=cols1[[2]], lwd=0.5, lty=3)
lines(xS, unlist(cts2_b_high[startingDay:length(cases1)]), col=cols1[[2]], lwd=0.5, lty=3)
slidingWindow = 7
for (i in (slidingWindow+1):length(cases1))
	{
		cases2[[i]] = mean(unlist(cases1[(i-slidingWindow):i]), na.rm=T)
		cts2_a[[i]] = -mean(unlist(cts1_a[(i-slidingWindow):i]), na.rm=T)
		cts2_a_low[[i]] = -mean(unlist(cts1_a_low[(i-slidingWindow):i]), na.rm=T)
		cts2_a_high[[i]] = -mean(unlist(cts1_a_high[(i-slidingWindow):i]), na.rm=T)
		cts2_b[[i]] = -mean(unlist(cts1_b[(i-slidingWindow):i]), na.rm=T)
		cts2_b_low[[i]] = -mean(unlist(cts1_b_low[(i-slidingWindow):i]), na.rm=T)
		cts2_b_high[[i]] = -mean(unlist(cts1_b_high[(i-slidingWindow):i]), na.rm=T)
	}
plot(xS, unlist(cases2[startingDay:length(cases1)]), col=NA, lwd=0.5, type="l", axes=F, ann=F, xlim=c(min(xS),max(xS)))
yy_l = c(rep(0,length(xS)),rev(unlist(cases2[startingDay:length(cases1)])))
xx_l = c(xS,rev(xS)); polygon(xx_l, yy_l, col=cols2[[1]], border=0)
axis(side=1, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.02, col.axis="gray30", mgp=c(0,0.03,0), at=labels1, labels=labels2)
axis(side=2, lwd.tick=0.2, cex.axis=0.7, lwd=0.2, tck=-0.02, col.axis="gray30", mgp=c(0,0.20,0), at=seq(0,200,50), pos=labels1[3]-3)
par(new=T); yMin = min(c(unlist(cts2_a_low),unlist(cts2_a_high)), na.rm=T); yMax = max(c(unlist(cts2_a_low),unlist(cts2_a_high)), na.rm=T)
plot(xS, unlist(cts2_a[startingDay:length(cases1)]), col=cols1[[2]], lwd=1.5, type="l", axes=F, ann=F, xlim=c(min(xS),max(xS)), ylim=c(-30,-8))
yy_l = c(rep(unlist(cts2_a_high[startingDay:length(cases1)])),rev(unlist(cts2_a_low[startingDay:length(cases1)])))
xx_l = c(xS,rev(xS)); polygon(xx_l, yy_l, col=cols2[[2]], border=0)
axis(side=4, lwd.tick=0.2, cex.axis=0.7, lwd=0.2, tck=-0.02, col.axis=cols1[[1]], mgp=c(0,0.10,0), at=seq(-40,0,5), labels=seq(40,0,-5))
lines(xS, unlist(cts2_b[startingDay:length(cases1)]), col=cols1[[2]], lwd=1.5, lty=2)
lines(xS, unlist(cts2_b_low[startingDay:length(cases1)]), col=cols1[[2]], lwd=0.5, lty=3)
lines(xS, unlist(cts2_b_high[startingDay:length(cases1)]), col=cols1[[2]], lwd=0.5, lty=3)
dev.off()

dev.new(width=2.3, height=2); par(oma=c(0,0,0,0), mar=c(1.0,1.0,0.2,0), col="gray30", lwd=0.2); slidingWindow = 14
for (i in (slidingWindow+1):length(cases1))
	{
		cases2[[i]] = mean(unlist(cases1[(i-slidingWindow):i]), na.rm=T)
		cts2_a[[i]] = -mean(unlist(cts1_a[(i-slidingWindow):i]), na.rm=T)
		cts2_a_low[[i]] = -mean(unlist(cts1_a_low[(i-slidingWindow):i]), na.rm=T)
		cts2_a_high[[i]] = -mean(unlist(cts1_a_high[(i-slidingWindow):i]), na.rm=T)
		cts2_b[[i]] = -mean(unlist(cts1_b[(i-slidingWindow):i]), na.rm=T)
		cts2_b_low[[i]] = -mean(unlist(cts1_b_low[(i-slidingWindow):i]), na.rm=T)
		cts2_b_high[[i]] = -mean(unlist(cts1_b_high[(i-slidingWindow):i]), na.rm=T)
	}
mat = cbind(unlist(cases2), unlist(cts2_a), unlist(cts2_b)); gaps = seq(0,30)
cors_a = rep(NA, length(gaps)); cors_b = rep(NA, length(gaps))
for (i in 1:length(gaps))
	{
		indices1 = startingDay:(dim(mat)[1]-gaps[i]); indices2 = (startingDay+gaps[i]):dim(mat)[1]
		tmp = cbind(mat[indices1,2], mat[indices2,1]); tmp = tmp[which((!is.na(tmp[,1]))&(!is.na(tmp[,2]))),]
		cors_a[i] = cor(tmp[,2], tmp[,1], method="spearman")
		tmp = cbind(mat[indices1,3], mat[indices2,1]); tmp = tmp[which((!is.na(tmp[,1]))&(!is.na(tmp[,2]))),]
		cors_b[i] = cor(tmp[,2], tmp[,1], method="spearman") 
	}
plot(gaps, cors_b, axes=F, ann=F, type="l", col="gray30", lwd=0.7, lty=2); # points(gaps, cors_b, col="gray30", cex=0.5, pch=16)
lines(gaps, cors_a, col="gray30", lwd=0.7); points(gaps, cors_a, col="gray30", cex=0.5, pch=16)
axis(side=1, lwd.tick=0.2, cex.axis=0.7, lwd=0.2, tck=-0.03, col.axis="gray30", mgp=c(0,0.05,0), at=seq(0,35,5))
axis(side=2, lwd.tick=0.2, cex.axis=0.7, lwd=0.2, tck=-0.03, col.axis="gray30", mgp=c(0,0.25,0), at=seq(0.45,0.85,0.10))

dev.new(width=2.3, height=2); par(oma=c(0,0,0,0), mar=c(1.0,1.0,0.2,0), col="gray30", lwd=0.2); slidingWindow = 7
for (i in (slidingWindow+1):length(cases1))
	{
		cases2[[i]] = mean(unlist(cases1[(i-slidingWindow):i]), na.rm=T)
		cts2_a[[i]] = -mean(unlist(cts1_a[(i-slidingWindow):i]), na.rm=T)
		cts2_a_low[[i]] = -mean(unlist(cts1_a_low[(i-slidingWindow):i]), na.rm=T)
		cts2_a_high[[i]] = -mean(unlist(cts1_a_high[(i-slidingWindow):i]), na.rm=T)
		cts2_b[[i]] = -mean(unlist(cts1_b[(i-slidingWindow):i]), na.rm=T)
		cts2_b_low[[i]] = -mean(unlist(cts1_b_low[(i-slidingWindow):i]), na.rm=T)
		cts2_b_high[[i]] = -mean(unlist(cts1_b_high[(i-slidingWindow):i]), na.rm=T)
	}
mat = cbind(unlist(cases2), unlist(cts2_a), unlist(cts2_b)); gaps = seq(0,30)
cors_a = rep(NA, length(gaps)); cors_b = rep(NA, length(gaps))
for (i in 1:length(gaps))
	{
		indices1 = startingDay:(dim(mat)[1]-gaps[i]); indices2 = (startingDay+gaps[i]):dim(mat)[1]
		tmp = cbind(mat[indices1,2], mat[indices2,1]); tmp = tmp[which((!is.na(tmp[,1]))&(!is.na(tmp[,2]))),]
		cors_a[i] = cor(tmp[,2], tmp[,1], method="spearman")
		tmp = cbind(mat[indices1,3], mat[indices2,1]); tmp = tmp[which((!is.na(tmp[,1]))&(!is.na(tmp[,2]))),]
		cors_b[i] = cor(tmp[,2], tmp[,1], method="spearman") 
	}
plot(gaps, cors_b, axes=F, ann=F, type="l", col="gray30", lwd=0.7, lty=2); # points(gaps, cors_b, col="gray30", cex=0.5, pch=16)
lines(gaps, cors_a, col="gray30", lwd=0.7); points(gaps, cors_a, col="gray30", cex=0.5, pch=16)
axis(side=1, lwd.tick=0.2, cex.axis=0.7, lwd=0.2, tck=-0.03, col.axis="gray30", mgp=c(0,0.05,0), at=seq(0,35,5))
axis(side=2, lwd.tick=0.2, cex.axis=0.7, lwd=0.2, tck=-0.03, col.axis="gray30", mgp=c(0,0.25,0), at=seq(0.45,0.85,0.10))

