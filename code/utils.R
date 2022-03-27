library(tidyverse) 
library(purrr)
library(ggrepel)

getslopes <- function(props_tidy, whichdrug="GSK", sgRNA_ID_special="EZH2_030", whichweekmap=weekmap){

	# Store the proportions for these special guides: 
	special_join_df <- props_tidy %>% 
			filter(sgRNA_ID==sgRNA_ID_special & drug==whichdrug) %>% 
			rename(sgRNA_ID_special=sgRNA_ID, prop_special=prop)
		
	# Calculate the difference in Z slope versus reference for each guide: 
	slope_df <- props_tidy %>% 
		left_join(whichweekmap,by="event") %>% 
		# Join proportions of "special" guides:
		inner_join(special_join_df, by=c("drug","event","rep")) %>% 
		# Calculate "Z", the thing whose slope gives the on/off drug growth rate:
		mutate(Z = log(prop)-log(prop_special)) %>% 
		# Format data frame for calculating slopes: 
		select(sgRNA_ID, drug, rep, event, week, Z) %>% 
		arrange(sgRNA_ID, drug, rep, week) %>% 
		group_by(sgRNA_ID, drug, rep) %>% 
		# Calculate slope of Z line at weeks 5 and 8: 
		mutate(slope=(Z-lag(Z))/(week-lag(week))) %>% 
		# Ged rid of Week 0:
		filter(!is.na(slope)) %>% 
		select(sgRNA_ID, drug, event, rep, slope)

	return(slope_df)
}

getslopes_uncertainty <- function(props_tidy, whichdrug="GSK", sgRNA_ID_special="EZH2_030", whichweekmap=weekmap, nsamp=10000, ndraw=1000){

	# Store the proportions for these special guides: 
	special_join_df <- props_tidy %>% 
			filter(sgRNA_ID==sgRNA_ID_special & drug==whichdrug) %>% 
			mutate(rowid=1:n()) %>% 
			split(.$rowid) %>% 
			map(~ tibble(sgRNA_ID=.$sgRNA_ID, drug=.$drug, event=.$event, rep=.$rep, prop=rbeta(ndraw,.$prop*nsamp+1, nsamp-.$prop*nsamp+1), draw=1:ndraw)) %>% 
			bind_rows() %>% 
			rename(sgRNA_ID_special=sgRNA_ID, prop_special=prop)
		
	# Calculate the difference in Z slope versus reference for each guide: 
	slope_df <- props_tidy %>% 
		left_join(whichweekmap,by="event") %>% 
		# Join proportions of "special" guides:
		inner_join(special_join_df, by=c("drug","event","rep","draw")) %>% 
		# Calculate "Z", the thing whose slope gives the on/off drug growth rate:
		mutate(Z = log(prop)-log(prop_special)) %>% 
		# Format data frame for calculating slopes: 
		select(draw, sgRNA_ID, drug, rep, event, week, Z) %>% 
		arrange(draw, sgRNA_ID, drug, rep, week) %>% 
		group_by(draw, sgRNA_ID, drug, rep) %>% 
		# Calculate slope of Z line at weeks 5 and 8: 
		mutate(slope=(Z-lag(Z))/(week-lag(week))) %>% 
		# Ged rid of Week 0:
		filter(!is.na(slope)) %>% 
		select(draw, sgRNA_ID, drug, event, rep, slope)

	return(slope_df)
}

# temp %>% 
# 	arrange(draw) %>% 
# 	select(sgRNA_ID, drug, event, rep, prop, draw, prop_special) %>% 
# 	print(n=50)


getrdiff <- function(props_tidy, whichdrug="GSK", sgRNA_ID_special="EZH2_030", whichweekmap=weekmap){

	slope_df <- getslopes(props_tidy, whichdrug, sgRNA_ID_special, whichweekmap)

	rdiff_df <- slope_df %>% 
		# Calculate change in slope from start-switch to switch-end:
		group_by(sgRNA_ID, drug, rep) %>% 
		mutate(rdiff=slope-lag(slope)) %>% 
		# Format data frame:
		filter(!is.na(rdiff)) %>% 
		select(sgRNA_ID, drug, rep, rdiff)

	return(rdiff_df)

}

getrdiff_uncertainty <- function(props_tidy, whichdrug="GSK", sgRNA_ID_special="EZH2_030", whichweekmap=weekmap, nsamp=10000, ndraw=1000){

	props_tidy_uncertainty <- props_tidy %>% 
		filter(drug==whichdrug) %>% 
		mutate(rowid=1:n()) %>% 
		split(.$rowid) %>% 
		map(~ tibble(sgRNA_ID=.$sgRNA_ID, drug=.$drug, event=.$event, rep=.$rep, prop_raw=.$prop, prop=rbeta(ndraw,.$prop*nsamp+1, nsamp-.$prop*nsamp+1), draw=1:ndraw)) %>% 
		bind_rows()

	# slope_df_uncertainty <- props_tidy_uncertainty %>% 
	# 	split(.$draw) %>% 
	# 	map(~ getslopes(., whichdrug, sgRNA_ID_special, whichweekmap)) %>% 
	# 	bind_rows(.id="draw") %>% 
	# 	mutate(draw=as.integer(draw))

	slope_df_uncertainty <- props_tidy_uncertainty %>% 
		split(.$draw) %>% 
		map(~ getslopes_uncertainty(., whichdrug, sgRNA_ID_special, whichweekmap, nsamp, ndraw)) %>% 
		bind_rows() %>% 
		mutate(draw=as.integer(draw))

	rdiff_df_uncertainty <- slope_df_uncertainty %>% 
		# Calculate change in slope from start-switch to switch-end:
		group_by(draw, sgRNA_ID, drug, rep) %>% 
		mutate(rdiff=slope-lag(slope)) %>% 
		# Format data frame:
		filter(!is.na(rdiff)) %>% 
		select(draw, sgRNA_ID, drug, rep, rdiff) %>% 
		ungroup() 

	return(rdiff_df_uncertainty)

}


summarise_rdiff_uncertainty <- function(rdiff_df_uncertainty){
	
	out <- rdiff_df_uncertainty %>% 
		# group_by(sgRNA_ID, drug) %>% 
		group_by(sgRNA_ID, drug, rep) %>% 
		summarise(
			rdiff_lwr95=quantile(rdiff,0.025),
			rdiff_lwr90=quantile(rdiff,0.05),
			rdiff_lwrQ=quantile(rdiff,0.25),
			rdiff_med=quantile(rdiff,0.5),
			rdiff_mean=mean(rdiff),
			rdiff_uprQ=quantile(rdiff,0.75),
			rdiff_upr90=quantile(rdiff,0.95),
			rdiff_upr95=quantile(rdiff,0.975)
			) %>% 
		ungroup() 

	return(out)
}

plot_growthcurves <- function(props_tidy, whichrep=0, whichdrug="GSK", sgRNA_ID_special="EZH2_030", bigval=0.1, bigspots=c("switch"), whichweekmap=weekmap, whichguideattributes=guide_attributes){

	slope_df <- getslopes(props_tidy, whichdrug, sgRNA_ID_special, whichweekmap)
	rdiff_df <- getrdiff(props_tidy, whichdrug, sgRNA_ID_special, whichweekmap)

	inputdf_on <- props_tidy %>% 
		filter(rep==whichrep & event=="start" & drug==whichdrug) %>% 
		select(sgRNA_ID, drug, propstart=prop) %>% 
		left_join(filter(slope_df, rep==whichrep & event=="switch" & drug==whichdrug), by=c("drug","sgRNA_ID"))

	inputdf_off <- props_tidy_onoff %>% 
		filter(rep==whichrep & event=="switch" & drug==whichdrug) %>% 
		select(sgRNA_ID, drug, propswitch=prop) %>% 
		left_join(filter(slope_df, rep==whichrep & event=="end" & drug==whichdrug), by=c("drug","sgRNA_ID"))

	rvec_on <- pull(inputdf_on, slope)
	rvec_off <- pull(inputdf_off, slope)
	x0vec_on <- pull(inputdf_on, propstart)
	x0vec_off <- pull(inputdf_off, propswitch)

	lkmod <- function(t,x,r){

		dx <- rep(0, length(x))
		for(i in 1:length(x)){
			dx[i] <- x[i]*(r[i] - sum(x*r))	
		}

		return(list(dx))

	}

	sol_on <- lsoda(y=x0vec_on, times=seq(from=0,to=5,by=0.01), func=lkmod, parms=rvec_on) %>% 
		as_tibble() %>% 
		pivot_longer(-time) %>% 
		mutate(name=as.numeric(name)) %>% 
		left_join((select(inputdf_on,sgRNA_ID) %>% mutate(name=1:n())),by="name") %>% 
		select(sgRNA_ID, time, value)

	sol_off <- lsoda(y=x0vec_off, times=seq(from=5,to=8,by=0.01), func=lkmod, parms=rvec_off) %>% 
		as_tibble() %>% 
		pivot_longer(-time) %>% 
		mutate(name=as.numeric(name)) %>% 
		left_join((select(inputdf_off,sgRNA_ID) %>% mutate(name=1:n())),by="name") %>% 
		select(sgRNA_ID, time, value)

	props_tidy_toplot <- left_join(props_tidy,whichweekmap,by="event") %>% 
		filter(rep==whichrep & drug==whichdrug) 

	bigguides <- props_tidy_toplot %>% 
		filter(event%in%bigspots & prop>bigval) %>% 
		pull(sgRNA_ID) %>% 
		unique()

	bigguides_addicted <- rdiff_df %>% 
		filter(drug==whichdrug & rep==whichrep & rdiff<0 & sgRNA_ID%in%bigguides) %>% 
		pull(sgRNA_ID) %>% 
		unique()

	bigguides_addicted_labs <- whichguideattributes %>% 
		filter(sgRNA_ID%in%bigguides_addicted) %>% 
		select(sgRNA_ID, lab) %>% 
		inner_join(
			(props_tidy_toplot %>% 
			filter(sgRNA_ID %in% bigguides_addicted) %>% 
			filter(event%in%c("switch","end"))), by="sgRNA_ID") %>% 
		select(sgRNA_ID, lab, prop, week)

	bigguides_unaddicted <- rdiff_df %>% 
		filter(drug==whichdrug & rep==whichrep & rdiff>0 & sgRNA_ID%in%bigguides) %>% 
		pull(sgRNA_ID) %>% 
		unique()

	bigguides_unaddicted_labs <- whichguideattributes %>% 
		filter(sgRNA_ID%in%bigguides_unaddicted) %>% 
		select(sgRNA_ID, lab) %>% 
		inner_join(
			(props_tidy_toplot %>% 
			filter(sgRNA_ID %in% bigguides_unaddicted) %>% 
			filter(event%in%c("switch","end"))), by="sgRNA_ID") %>% 
		select(sgRNA_ID, lab, prop, week)

	specialguides_labs <- whichguideattributes %>% 
		filter(sgRNA_ID==sgRNA_ID_special) %>% 
		select(sgRNA_ID, lab) %>% 
		inner_join(
			(props_tidy_toplot %>% 
			filter(sgRNA_ID==sgRNA_ID_special) %>% 
			filter(event%in%c("switch","end"))), by="sgRNA_ID") %>% 
		select(sgRNA_ID, lab, prop, week)

	fig_props <- ggplot() + 
		# BACKGROUND
		annotate("rect", xmin=filter(weekmap, event=="start")$week, xmax=filter(weekmap, event=="switch")$week, ymin=-Inf, ymax=Inf, fill="darkgray", alpha=0.2) +
		# THIN LINES
		geom_line(data=sol_on, aes(x=time, y=value, group=sgRNA_ID), col="darkgray", size=0.25, alpha=0.5) + 
		geom_line(data=sol_off, aes(x=time, y=value, group=sgRNA_ID), col="darkgray", size=0.25, alpha=0.5) + 
		# THICK LINES, ADDICTED
		geom_line(data=filter(sol_on,sgRNA_ID%in%bigguides_addicted), aes(x=time, y=value, group=sgRNA_ID), col="red", size=0.5) + 
		geom_line(data=filter(sol_off,sgRNA_ID%in%bigguides_addicted), aes(x=time, y=value, group=sgRNA_ID), col="red", size=0.5) + 
		# THICK LINES, UNADDICTED
		geom_line(data=filter(sol_on,sgRNA_ID%in%bigguides_unaddicted), aes(x=time, y=value, group=sgRNA_ID), col="blue", size=0.5) + 
		geom_line(data=filter(sol_off,sgRNA_ID%in%bigguides_unaddicted), aes(x=time, y=value, group=sgRNA_ID), col="blue", size=0.5) + 
		# THICK LINES, CONTROL
		geom_line(data=filter(sol_on,sgRNA_ID==sgRNA_ID_special), aes(x=time, y=value, group=sgRNA_ID), col="black", size=0.5) + 
		geom_line(data=filter(sol_off,sgRNA_ID==sgRNA_ID_special), aes(x=time, y=value, group=sgRNA_ID), col="black", size=0.5) + 
		# THIN POINTS
		geom_point(data=filter(props_tidy_toplot, event=="start"),aes(x=week, y=prop), col="darkgray", size=0.5, alpha=0.5)  + 
		geom_point(data=filter(props_tidy_toplot, event=="switch"),aes(x=week,y=prop), col="darkgray", size=0.5, alpha=0.5)  + 
		geom_point(data=filter(props_tidy_toplot, event=="end"),aes(x=week,y=prop), col="darkgray", size=0.5, alpha=0.5)  + 
		# THICK POINTS, ADDICTED
		geom_point(data=filter(props_tidy_toplot, sgRNA_ID%in%bigguides_addicted & event=="start"),aes(x=week, y=prop), col="red", size=1)  + 
		geom_point(data=filter(props_tidy_toplot, sgRNA_ID%in%bigguides_addicted & event=="switch"),aes(x=week,y=prop), col="red", size=1)  + 
		geom_point(data=filter(props_tidy_toplot, sgRNA_ID%in%bigguides_addicted & event=="end"),aes(x=week,y=prop), col="red", size=1)  + 
		# THICK POINTS, UNADDICTED
		geom_point(data=filter(props_tidy_toplot, sgRNA_ID%in%bigguides_unaddicted & event=="start"),aes(x=week, y=prop), col="blue", size=1)  + 
		geom_point(data=filter(props_tidy_toplot, sgRNA_ID%in%bigguides_unaddicted & event=="switch"),aes(x=week,y=prop), col="blue", size=1)  + 
		geom_point(data=filter(props_tidy_toplot, sgRNA_ID%in%bigguides_unaddicted & event=="end"),aes(x=week,y=prop), col="blue", size=1)  +
		# THICK POINTS, CONTROL
		geom_point(data=filter(props_tidy_toplot, sgRNA_ID==sgRNA_ID_special & event=="start"),aes(x=week, y=prop), col="black", size=1)  + 
		geom_point(data=filter(props_tidy_toplot, sgRNA_ID==sgRNA_ID_special & event=="switch"),aes(x=week,y=prop), col="black", size=1)  + 
		geom_point(data=filter(props_tidy_toplot, sgRNA_ID==sgRNA_ID_special & event=="end"),aes(x=week,y=prop), col="black", size=1)  +
		# LABELS 
		# geom_text(data=bigguides_addicted_labs, aes(x=week, y=prop, label=lab),  col="red", size=2, hjust=-0.1, vjust=0) + 
		# geom_text(data=bigguides_unaddicted_labs, aes(x=week, y=prop, label=lab),  col="blue", size=2, hjust=-0.1, vjust=0) +
		# geom_text(data=specialguides_labs, aes(x=week, y=prop, label=lab),  col="black", size=2, hjust=-0.1, vjust=0) + 
		geom_text_repel(data=bigguides_addicted_labs, aes(x=week, y=prop, label=lab),  col="red", size=2, min.segment.length=0.1, segment.size=0.3, box.padding=0.2) + 
		geom_text_repel(data=bigguides_unaddicted_labs, aes(x=week, y=prop, label=lab),  col="blue", size=2, min.segment.length=0.1, segment.size=0.3, box.padding=0.2) +
		geom_text_repel(data=specialguides_labs, aes(x=week, y=prop, label=lab),  col="black", size=2, min.segment.length=0.1, segment.size=0.3, box.padding=0.2) + 
		# PLOT ATTRIBUTES 
		scale_x_continuous(limits=c(min(whichweekmap$week),max(whichweekmap$week)+0.5), breaks=min(whichweekmap$week):(max(whichweekmap$week)), minor_breaks=min(whichweekmap$week):(max(whichweekmap$week))) + 
		theme_classic() + 
		theme(legend.position="none") + 
		labs(title=paste0(whichdrug,": Replicate ",whichrep), x="Week", y="Proportion")

	return(fig_props)


}




plot_growthcurves_uncertainty <- function(props_tidy, rdiff_df_uncertainty_drug_summary, whichrep=0, whichdrug="GSK", sgRNA_ID_special="EZH2_030", bigval=0.1, bigspots=c("switch"), whichweekmap=weekmap, whichguideattributes=guide_attributes){

	slope_df <- getslopes(props_tidy, whichdrug, sgRNA_ID_special, whichweekmap)

	inputdf_on <- props_tidy %>% 
		filter(rep==whichrep & event=="start" & drug==whichdrug) %>% 
		select(sgRNA_ID, drug, propstart=prop) %>% 
		left_join(filter(slope_df, rep==whichrep & event=="switch" & drug==whichdrug), by=c("drug","sgRNA_ID"))

	inputdf_off <- props_tidy_onoff %>% 
		filter(rep==whichrep & event=="switch" & drug==whichdrug) %>% 
		select(sgRNA_ID, drug, propswitch=prop) %>% 
		left_join(filter(slope_df, rep==whichrep & event=="end" & drug==whichdrug), by=c("drug","sgRNA_ID"))

	rvec_on <- pull(inputdf_on, slope)
	rvec_off <- pull(inputdf_off, slope)
	x0vec_on <- pull(inputdf_on, propstart)
	x0vec_off <- pull(inputdf_off, propswitch)

	lkmod <- function(t,x,r){

		dx <- rep(0, length(x))
		for(i in 1:length(x)){
			dx[i] <- x[i]*(r[i] - sum(x*r))	
		}

		return(list(dx))

	}

	sol_on <- lsoda(y=x0vec_on, times=seq(from=0,to=5,by=0.01), func=lkmod, parms=rvec_on) %>% 
		as_tibble() %>% 
		pivot_longer(-time) %>% 
		mutate(name=as.numeric(name)) %>% 
		left_join((select(inputdf_on,sgRNA_ID) %>% mutate(name=1:n())),by="name") %>% 
		select(sgRNA_ID, time, value)

	sol_off <- lsoda(y=x0vec_off, times=seq(from=5,to=8,by=0.01), func=lkmod, parms=rvec_off) %>% 
		as_tibble() %>% 
		pivot_longer(-time) %>% 
		mutate(name=as.numeric(name)) %>% 
		left_join((select(inputdf_off,sgRNA_ID) %>% mutate(name=1:n())),by="name") %>% 
		select(sgRNA_ID, time, value)

	props_tidy_toplot <- left_join(props_tidy,whichweekmap,by="event") %>% 
		filter(rep==whichrep & drug==whichdrug) 

	bigguides <- props_tidy_toplot %>% 
		filter(event%in%bigspots & prop>bigval) %>% 
		pull(sgRNA_ID) %>% 
		unique()

	allguides_addicted <- rdiff_df_uncertainty_drug_summary %>% 
		filter(rep %in% whichrep) %>% 
		filter(rdiff_upr90<0) %>% 
		pull(sgRNA_ID) %>% 
		unique()

	bigguides_addicted <- intersect(bigguides, allguides_addicted) 

	allguides_addicted_labs <- whichguideattributes %>% 
		filter(sgRNA_ID%in%allguides_addicted) %>% 
		select(sgRNA_ID, lab) %>% 
		inner_join(
			(props_tidy_toplot %>% 
			filter(sgRNA_ID %in% allguides_addicted) %>% 
			filter(event%in%c("switch","end"))), by="sgRNA_ID") %>% 
		select(sgRNA_ID, lab, prop, week)

	bigguides_addicted_labs <- whichguideattributes %>% 
		filter(sgRNA_ID%in%bigguides_addicted) %>% 
		select(sgRNA_ID, lab) %>% 
		inner_join(
			(props_tidy_toplot %>% 
			filter(sgRNA_ID %in% bigguides_addicted) %>% 
			filter(event%in%c("switch","end"))), by="sgRNA_ID") %>% 
		select(sgRNA_ID, lab, prop, week)

	allguides_unaddicted <- rdiff_df_uncertainty_drug_summary %>% 
		filter(rep %in% whichrep) %>% 
		filter(rdiff_mean>0) %>% 
		pull(sgRNA_ID) %>% 
		unique()

	bigguides_unaddicted <- intersect(bigguides, allguides_unaddicted) 

	bigguides_unaddicted_labs <- whichguideattributes %>% 
		filter(sgRNA_ID%in%bigguides_unaddicted) %>% 
		select(sgRNA_ID, lab) %>% 
		inner_join(
			(props_tidy_toplot %>% 
			filter(sgRNA_ID %in% bigguides_unaddicted) %>% 
			filter(event%in%c("switch","end"))), by="sgRNA_ID") %>% 
		select(sgRNA_ID, lab, prop, week)

	specialguides_labs <- whichguideattributes %>% 
		filter(sgRNA_ID==sgRNA_ID_special) %>% 
		select(sgRNA_ID, lab) %>% 
		inner_join(
			(props_tidy_toplot %>% 
			filter(sgRNA_ID==sgRNA_ID_special) %>% 
			filter(event%in%c("switch","end"))), by="sgRNA_ID") %>% 
		select(sgRNA_ID, lab, prop, week)

	fig_props <- ggplot() + 
		# BACKGROUND
		annotate("rect", xmin=filter(weekmap, event=="start")$week, xmax=filter(weekmap, event=="switch")$week, ymin=-Inf, ymax=Inf, fill="darkgray", alpha=0.2) +
		# THIN LINES
		geom_line(data=sol_on, aes(x=time, y=value, group=sgRNA_ID), col="darkgray", size=0.25, alpha=0.5) + 
		geom_line(data=sol_off, aes(x=time, y=value, group=sgRNA_ID), col="darkgray", size=0.25, alpha=0.5) + 
		# THICK LINES, BIG ADDICTED
		geom_line(data=filter(sol_on,sgRNA_ID%in%bigguides_addicted), aes(x=time, y=value, group=sgRNA_ID), col="red", size=0.5) + 
		geom_line(data=filter(sol_off,sgRNA_ID%in%bigguides_addicted), aes(x=time, y=value, group=sgRNA_ID), col="red", size=0.5) + 
		# THICK LINES, BIG UNADDICTED
		geom_line(data=filter(sol_on,sgRNA_ID%in%bigguides_unaddicted), aes(x=time, y=value, group=sgRNA_ID), col="blue", size=0.5) + 
		geom_line(data=filter(sol_off,sgRNA_ID%in%bigguides_unaddicted), aes(x=time, y=value, group=sgRNA_ID), col="blue", size=0.5) + 
		# THICK LINES, CONTROL
		geom_line(data=filter(sol_on,sgRNA_ID==sgRNA_ID_special), aes(x=time, y=value, group=sgRNA_ID), col="black", size=0.5) + 
		geom_line(data=filter(sol_off,sgRNA_ID==sgRNA_ID_special), aes(x=time, y=value, group=sgRNA_ID), col="black", size=0.5) + 
		# THIN POINTS
		geom_point(data=filter(props_tidy_toplot, event=="start"),aes(x=week, y=prop), col="darkgray", size=0.5, alpha=0.5)  + 
		geom_point(data=filter(props_tidy_toplot, event=="switch"),aes(x=week,y=prop), col="darkgray", size=0.5, alpha=0.5)  + 
		geom_point(data=filter(props_tidy_toplot, event=="end"),aes(x=week,y=prop), col="darkgray", size=0.5, alpha=0.5)  + 
		# THICK POINTS, BIG ADDICTED
		geom_point(data=filter(props_tidy_toplot, sgRNA_ID%in%bigguides_addicted & event=="start"),aes(x=week, y=prop), col="red", size=1)  + 
		geom_point(data=filter(props_tidy_toplot, sgRNA_ID%in%bigguides_addicted & event=="switch"),aes(x=week,y=prop), col="red", size=1)  + 
		geom_point(data=filter(props_tidy_toplot, sgRNA_ID%in%bigguides_addicted & event=="end"),aes(x=week,y=prop), col="red", size=1)  + 
		# THICK POINTS, BIG UNADDICTED
		geom_point(data=filter(props_tidy_toplot, sgRNA_ID%in%bigguides_unaddicted & event=="start"),aes(x=week, y=prop), col="blue", size=1)  + 
		geom_point(data=filter(props_tidy_toplot, sgRNA_ID%in%bigguides_unaddicted & event=="switch"),aes(x=week,y=prop), col="blue", size=1)  + 
		geom_point(data=filter(props_tidy_toplot, sgRNA_ID%in%bigguides_unaddicted & event=="end"),aes(x=week,y=prop), col="blue", size=1)  +
		# THICK POINTS, CONTROL
		geom_point(data=filter(props_tidy_toplot, sgRNA_ID==sgRNA_ID_special & event=="start"),aes(x=week, y=prop), col="black", size=1)  + 
		geom_point(data=filter(props_tidy_toplot, sgRNA_ID==sgRNA_ID_special & event=="switch"),aes(x=week,y=prop), col="black", size=1)  + 
		geom_point(data=filter(props_tidy_toplot, sgRNA_ID==sgRNA_ID_special & event=="end"),aes(x=week,y=prop), col="black", size=1)  +
		# LABELS 
		# geom_text(data=bigguides_addicted_labs, aes(x=week, y=prop, label=lab),  col="red", size=2, hjust=-0.1, vjust=0) + 
		# geom_text(data=bigguides_unaddicted_labs, aes(x=week, y=prop, label=lab),  col="blue", size=2, hjust=-0.1, vjust=0) +
		# geom_text(data=specialguides_labs, aes(x=week, y=prop, label=lab),  col="black", size=2, hjust=-0.1, vjust=0) + 
		geom_text_repel(data=bigguides_addicted_labs, aes(x=week, y=prop, label=lab),  col="red", size=2, min.segment.length=0.1, segment.size=0.3, box.padding=0.2) + 
		geom_text_repel(data=bigguides_unaddicted_labs, aes(x=week, y=prop, label=lab),  col="blue", size=2, min.segment.length=0.1, segment.size=0.3, box.padding=0.2) +
		geom_text_repel(data=specialguides_labs, aes(x=week, y=prop, label=lab),  col="black", size=2, min.segment.length=0.1, segment.size=0.3, box.padding=0.2) + 
		# PLOT ATTRIBUTES 
		scale_x_continuous(limits=c(min(whichweekmap$week),max(whichweekmap$week)+0.5), breaks=min(whichweekmap$week):(max(whichweekmap$week)), minor_breaks=min(whichweekmap$week):(max(whichweekmap$week))) + 
		theme_classic() + 
		theme(legend.position="none") + 
		labs(title=paste0(whichdrug,": Replicate ",whichrep), x="Week", y="Proportion")

	return(fig_props)


}

