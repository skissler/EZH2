# ==============================================================================
# Import
# ==============================================================================

print(paste0("Starting script: ",Sys.time()))

# Load packages:
library(tidyverse) 
library(ggrepel)
library(purrr)
library(deSolve)
source('code/utils.R') 

# Set key parameters: 
bigpropcutoff <- 0.05 # Minimum poportion for addiction cutoff
weekmap <- tibble(week=c(0,5,8),event=c("start","switch","end")) %>% 
	mutate(event=factor(event, levels=c("start","switch","end")))
special_guide_GSK <- "EZH2_030"
special_guide_EED <- "EED_099"

# Import raw proportions data: 
props_raw <- read_csv("data/AddictionProp_031022.csv",
	col_types=list(  
		col_character(),		# "sgRNA_ID"
		col_character(),		# "Gene"
		col_character(),		# "sgRNA_seq"
		col_double(),			# "mitSpecScore"
		col_double(),			# "DoenchScore"
		col_double(),			# "cut_site_NT"
		col_double(),			# "cut_site_AA"
		col_character(),		# "Domain"
		col_double(),			# "d0"
		col_double(),			# "w5_GSKR1"
		col_double(),			# "w5_GSKR2"
		col_double(),			# "w5_GSKR3"
		col_double(),			# "w5_EEDR1"
		col_double(),			# "w5_EEDR2"
		col_double(),			# "w5_EEDR3"
		col_double(),			# "w8_GSKR1off"
		col_double(),			# "w8_GSKR2off"
		col_double(),			# "w8_GSKR3off"
		col_double(),			# "w8_EEDR1off"
		col_double(),			# "w8_EEDR2off"
		col_double(),			# "w8_EEDR3off"
		col_double(),			# "w8_GSKR1on"
		col_double(),			# "w8_GSKR2on"
		col_double(),			# "w8_GSKR3on"
		col_double(),			# "w8_EEDR1on"
		col_double(),			# "w8_EEDR2on"
		col_double()),			# "w8_EEDR3on"
	col_names=c(
		"sgRNA_ID",
		"Gene",
		"sgRNA_seq",
		"mitSpecScore",
		"DoenchScore",
		"cut_site_NT",
		"cut_site_AA",
		"Domain",
		"d0",
		"w5_GSKR1",
		"w5_GSKR2",
		"w5_GSKR3",
		"w5_EEDR1",
		"w5_EEDR2",
		"w5_EEDR3",
		"w8_GSKR1off",
		"w8_GSKR2off",
		"w8_GSKR3off",
		"w8_EEDR1off",
		"w8_EEDR2off",
		"w8_EEDR3off",
		"w8_GSKR1on",
		"w8_GSKR2on",
		"w8_GSKR3on",
		"w8_EEDR1on",
		"w8_EEDR2on",
		"w8_EEDR3on"),
	skip=1)

# Store guide attributes: 
guide_attributes <- props_raw %>% 
	select(sgRNA_ID, Gene, sgRNA_seq, mitSpecScore, DoenchScore, cut_site_NT, cut_site_AA, Domain) %>% 
	mutate(lab=paste0(Gene," ",cut_site_AA))


# Put guide proportions into tidy format: 
# Start by formatting the day-0 data: 
props_d0 <- props_raw %>% 
	select(sgRNA_ID, prop=d0) %>% 
	split(.$sgRNA_ID) %>% 
	# repeat day-0 proportions across drugs and reps: 
	map(~ bind_rows(
		mutate(., drug="GSK", condition="on", week=0, rep=1),
		mutate(., drug="GSK", condition="on", week=0, rep=2),
		mutate(., drug="GSK", condition="on", week=0, rep=3),
		mutate(., drug="EED", condition="on", week=0, rep=1),
		mutate(., drug="EED", condition="on", week=0, rep=2),
		mutate(., drug="EED", condition="on", week=0, rep=3)
		))  %>% 
	bind_rows() %>% 
	select(sgRNA_ID, drug, condition, week, rep, prop)

# Now format the post-day-0 data: 
props_postd0 <- props_raw %>% 
	# Keep only the sgRNA_ID and the proportions: 
	select(
		sgRNA_ID, 
		w5_GSKR1,
		w5_GSKR2,
		w5_GSKR3,
		w5_EEDR1,
		w5_EEDR2,
		w5_EEDR3,
		w8_GSKR1off,
		w8_GSKR2off,
		w8_GSKR3off,
		w8_EEDR1off,
		w8_EEDR2off,
		w8_EEDR3off,
		w8_GSKR1on,
		w8_GSKR2on,
		w8_GSKR3on,
		w8_EEDR1on,
		w8_EEDR2on,
		w8_EEDR3on) %>% 
	# Put in long format: 
	pivot_longer(-sgRNA_ID, values_to="prop") %>% 
	# Mark whether the timepoint is at the end of an on-drug or off-drug period: 
	mutate(condition=case_when(
		substr(name,9,1000) %in% c("","on")~"on",TRUE~"off")) %>% 
	# Store the drug name: 
	mutate(drug=substr(name,4,6)) %>% 
	mutate(drug=case_when(drug==""~NA_character_,TRUE~drug)) %>% 
	# Store the replicate number: 
	mutate(rep=as.numeric(substr(name,8,8))) %>% 
	# Store the week number:
	mutate(week=as.numeric(substr(name,2,2))) %>% 
	# Toss the original field names:
	select(sgRNA_ID, drug, condition, week, rep, prop) 

# Now combine day-0 and post-day-0 data into a tidy data frame: 
props_tidy_full <- bind_rows(props_d0, props_postd0) %>% 
	# Add 1e-12 to all proportions to avoid 0s:
	mutate(prop=prop+1e-12) %>%
	group_by(drug, condition, week, rep) %>% 
	# Re-normalize proportions to ensure they sum to 1:
	mutate(prop=prop/sum(prop)) %>% 
	ungroup() %>% 
	left_join(select(guide_attributes, sgRNA_ID, Gene, cut_site_AA), by="sgRNA_ID") %>% 
	mutate(cut_site_AA=case_when(cut_site_AA<0~1e12, TRUE~cut_site_AA)) %>% 
	arrange(drug, Gene, cut_site_AA, week, desc(condition), rep) %>% 
	select(-Gene, -cut_site_AA)

# Append average proportions marked with rep=0: 
props_tidy_avg <- props_tidy_full %>% 
	group_by(sgRNA_ID, drug, condition, week) %>% 
	summarise(prop=mean(prop)) %>% 
	ungroup() %>% 
	mutate(rep=0) %>% 
	select(sgRNA_ID, drug, condition, week, rep, prop)

# Add the average proportions onto props_tidy_full: 
props_tidy_full <- bind_rows(props_tidy_full, props_tidy_avg)
rm(props_tidy_avg)
# rm(props_tidy_avg_2)

# Add an 'events' column onto props_tidy_full: 
props_tidy_full <- props_tidy_full %>% 
	left_join(weekmap, by="week") %>% 
	select(sgRNA_ID, drug, condition, event, rep, prop) 
	
# Separate into data frames for the on/off experiment and the on/on experiment:
props_tidy_onoff <- props_tidy_full %>% 
	mutate(tocut=case_when(event=="end" & condition=="on"~1, TRUE~0)) %>% 
	filter(tocut==0) %>% 
	select(-tocut) %>% 
	select(-condition)

props_tidy_onon <- props_tidy_full %>% 
	mutate(tocut=case_when(event=="end" & condition=="off"~1, TRUE~0)) %>% 
	filter(tocut==0) %>% 
	select(-tocut) %>% 
	select(-condition)

print(paste0("Done with import: ",Sys.time()))

# ==============================================================================
# Plot growth curves
# ==============================================================================

fig_growthcurve_avg_GSK <- plot_growthcurves(props_tidy_onoff, whichrep=0, whichdrug="GSK", sgRNA_ID_special=special_guide_GSK, bigval=bigpropcutoff, bigspots=c("switch"), whichweekmap=weekmap, whichguideattributes=guide_attributes) + 
	labs(title="GSK") + 
	scale_x_continuous(limits=c(0,8), expand=c(0,0.001), breaks=0:8) + 
	scale_y_continuous(limits=c(0,0.6), expand=c(0,0.001), breaks=seq(from=0,to=1,by=0.2)) 
fig_growthcurve_R1_GSK <- plot_growthcurves(props_tidy_onoff, whichrep=1, whichdrug="GSK", sgRNA_ID_special=special_guide_GSK, bigval=bigpropcutoff, bigspots=c("switch"), whichweekmap=weekmap, whichguideattributes=guide_attributes)
fig_growthcurve_R2_GSK <- plot_growthcurves(props_tidy_onoff, whichrep=2, whichdrug="GSK", sgRNA_ID_special=special_guide_GSK, bigval=bigpropcutoff, bigspots=c("switch"), whichweekmap=weekmap, whichguideattributes=guide_attributes)
fig_growthcurve_R3_GSK <- plot_growthcurves(props_tidy_onoff, whichrep=3, whichdrug="GSK", sgRNA_ID_special=special_guide_GSK, bigval=bigpropcutoff, bigspots=c("switch"), whichweekmap=weekmap, whichguideattributes=guide_attributes)

ggsave(fig_growthcurve_avg_GSK,file="figures/growthcurve_avg_GSK.pdf", width=2.5, height=2.5)
ggsave(fig_growthcurve_avg_GSK,file="figures/growthcurve_avg_GSK_rect.pdf", width=4, height=2.5)
ggsave(fig_growthcurve_avg_GSK,file="figures/growthcurve_avg_GSK.png", width=2.5, height=2.5)
ggsave(fig_growthcurve_avg_GSK,file="figures/growthcurve_avg_GSK_rect.png", width=4, height=2.5)
ggsave(fig_growthcurve_R1_GSK,file="figures/growthcurve_R1_GSK.pdf", width=2.5, height=2.5)
ggsave(fig_growthcurve_R1_GSK,file="figures/growthcurve_R1_GSK.png", width=2.5, height=2.5)
ggsave(fig_growthcurve_R2_GSK,file="figures/growthcurve_R2_GSK.pdf", width=2.5, height=2.5)
ggsave(fig_growthcurve_R2_GSK,file="figures/growthcurve_R2_GSK.png", width=2.5, height=2.5)
ggsave(fig_growthcurve_R3_GSK,file="figures/growthcurve_R3_GSK.pdf", width=2.5, height=2.5)
ggsave(fig_growthcurve_R3_GSK,file="figures/growthcurve_R3_GSK.png", width=2.5, height=2.5)


fig_growthcurve_avg_EED <- plot_growthcurves(props_tidy_onoff, whichrep=0, whichdrug="EED", sgRNA_ID_special=special_guide_EED, bigval=bigpropcutoff, bigspots=c("switch"), whichweekmap=weekmap, whichguideattributes=guide_attributes) + 
	labs(title="EED") + 
	scale_x_continuous(limits=c(0,8), expand=c(0,0.001), breaks=0:8) + 
	scale_y_continuous(limits=c(0,0.8), expand=c(0,0.001), breaks=seq(from=0,to=1,by=0.2)) 
fig_growthcurve_R1_EED <- plot_growthcurves(props_tidy_onoff, whichrep=1, whichdrug="EED", sgRNA_ID_special=special_guide_EED, bigval=bigpropcutoff, bigspots=c("switch"), whichweekmap=weekmap, whichguideattributes=guide_attributes)
fig_growthcurve_R2_EED <- plot_growthcurves(props_tidy_onoff, whichrep=2, whichdrug="EED", sgRNA_ID_special=special_guide_EED, bigval=bigpropcutoff, bigspots=c("switch"), whichweekmap=weekmap, whichguideattributes=guide_attributes)
fig_growthcurve_R3_EED <- plot_growthcurves(props_tidy_onoff, whichrep=3, whichdrug="EED", sgRNA_ID_special=special_guide_EED, bigval=bigpropcutoff, bigspots=c("switch"), whichweekmap=weekmap, whichguideattributes=guide_attributes)


ggsave(fig_growthcurve_avg_EED,file="figures/growthcurve_avg_EED.pdf", width=2.5, height=2.5)
ggsave(fig_growthcurve_avg_EED,file="figures/growthcurve_avg_EED_rect.pdf", width=4, height=2.5)
ggsave(fig_growthcurve_avg_EED,file="figures/growthcurve_avg_EED.png", width=2.5, height=2.5)
ggsave(fig_growthcurve_avg_EED,file="figures/growthcurve_avg_EED_rect.png", width=4, height=2.5)
ggsave(fig_growthcurve_R1_EED,file="figures/growthcurve_R1_EED.pdf", width=2.5, height=2.5)
ggsave(fig_growthcurve_R1_EED,file="figures/growthcurve_R1_EED.png", width=2.5, height=2.5)
ggsave(fig_growthcurve_R2_EED,file="figures/growthcurve_R2_EED.pdf", width=2.5, height=2.5)
ggsave(fig_growthcurve_R2_EED,file="figures/growthcurve_R2_EED.png", width=2.5, height=2.5)
ggsave(fig_growthcurve_R3_EED,file="figures/growthcurve_R3_EED.pdf", width=2.5, height=2.5)
ggsave(fig_growthcurve_R3_EED,file="figures/growthcurve_R3_EED.png", width=2.5, height=2.5)

print(paste0("Done with growth curves: ",Sys.time()))

# ==============================================================================
# Calculate rdiff
# ==============================================================================

rdiff_df_GSK <- getrdiff(props_tidy_onoff, whichdrug="GSK", sgRNA_ID_special=special_guide_GSK, whichweekmap=weekmap) %>% 
	# Calculate the average rdiff across replicates (late rdiff averaging):
	mutate(dummy=1) %>% 
	split(.$dummy) %>% 
	map(~ list(a=., b=(filter(., rep>0) %>% group_by(sgRNA_ID, drug) %>% summarise(rdiff=mean(rdiff), dummy=first(dummy)) %>% mutate(rep=-1) %>% select(sgRNA_ID, drug, rep, rdiff, dummy)))) %>% 
	map(~ bind_rows(.)) %>% 
	bind_rows() %>% 
	select(-dummy) %>% 
	arrange(sgRNA_ID, drug, rep)

# rdiff_df_uncertainty_GSK <- getrdiff_uncertainty(props_tidy_onoff, whichdrug="GSK", sgRNA_ID_special=special_guide_GSK, whichweekmap=weekmap, nsamp=10000, ndraw=1000)

print(paste0("Done with GSK rdiff: ",Sys.time()))

rdiff_df_EED <- getrdiff(props_tidy_onoff, whichdrug="EED", sgRNA_ID_special=special_guide_EED, whichweekmap=weekmap) %>% 
	# Calculate the average rdiff across replicates (late rdiff averaging):
	mutate(dummy=1) %>% 
	split(.$dummy) %>% 
	map(~ list(a=., b=(filter(., rep>0) %>% group_by(sgRNA_ID, drug) %>% summarise(rdiff=mean(rdiff), dummy=first(dummy)) %>% mutate(rep=-1) %>% select(sgRNA_ID, drug, rep, rdiff, dummy)))) %>% 
	map(~ bind_rows(.)) %>% 
	bind_rows() %>% 
	select(-dummy) %>% 
	arrange(sgRNA_ID, drug, rep)

# rdiff_df_uncertainty_EED <- getrdiff_uncertainty(props_tidy_onoff, whichdrug="EED", sgRNA_ID_special=special_guide_EED, whichweekmap=weekmap, nsamp=10000, ndraw=1000)

print(paste0("Done with EED rdiff: ",Sys.time()))

# ==============================================================================
# Identify addicted guides
# ==============================================================================

# Make a list of replicates with good control representation: 
validreps <- props_tidy_onoff %>% 
	filter((sgRNA_ID==special_guide_GSK & drug=="GSK") | (sgRNA_ID==special_guide_EED & drug=="EED")) %>% 
	# Rule: require ref guide to be >5% prevalence by "switch" time: 
	filter(event=="switch") %>% 
	filter(rep>0) %>% 
	filter(prop>bigpropcutoff) %>% 
	select(drug, rep) %>% 
	split(.$drug) %>% 
	map(~ .$rep) 

# Make a data frame of guides >5% at switch for each drug/replicate: 
bigswitchguides <- props_tidy_onoff %>% 
	# filter(rep>0) %>% 
	filter(event=="switch") %>% 
	filter(prop>bigpropcutoff) %>% 
	select(drug, rep, sgRNA_ID)

# Identify addicted guides for GSK: 
# rdiff_df_uncertainty_GSK_summary <- summarise_rdiff_uncertainty(rdiff_df_uncertainty_GSK)

addicted_df_GSK <- rdiff_df_GSK %>% 
	# filter(rep>0) %>% 
	inner_join(bigswitchguides, by=c("sgRNA_ID","drug","rep")) %>% 
	filter(rdiff<0) %>% 
	select(sgRNA_ID, drug, rep, rdiff) %>% 
	left_join(select(guide_attributes, sgRNA_ID, lab), by="sgRNA_ID")

# addicted_df_GSK_uncertainty <- rdiff_df_uncertainty_GSK_summary %>% 
# 	# filter(rep %in% validreps$GSK) %>% 
# 	filter(rep>0) %>% 
# 	inner_join(bigswitchguides, by=c("sgRNA_ID","drug","rep")) %>% 
# 	filter(rdiff_upr90<0) %>% 
# 	select(sgRNA_ID, drug, rep, rdiff_lwr90, rdiff_mean, rdiff_upr90) %>% 
# 	left_join(select(guide_attributes, sgRNA_ID, lab), by="sgRNA_ID")


# Identify addicted guides for EED: 
# rdiff_df_uncertainty_EED_summary <- summarise_rdiff_uncertainty(rdiff_df_uncertainty_EED)

addicted_df_EED <- rdiff_df_EED %>% 
	# filter(rep>0) %>% 
	inner_join(bigswitchguides, by=c("sgRNA_ID","drug","rep")) %>% 
	filter(rdiff<0) %>% 
	select(sgRNA_ID, drug, rep, rdiff) %>% 
	left_join(select(guide_attributes, sgRNA_ID, lab), by="sgRNA_ID")

# addicted_df_EED_uncertainty <- rdiff_df_uncertainty_EED_summary %>% 
# 	# filter(rep %in% validreps$EED) %>% 
# 	filter(rep>0) %>% 
# 	inner_join(bigswitchguides, by=c("sgRNA_ID","drug","rep")) %>% 
# 	filter(rdiff_upr90<0) %>% 
# 	select(sgRNA_ID, drug, rep, rdiff_lwr90, rdiff_mean, rdiff_upr90) %>% 
# 	left_join(select(guide_attributes, sgRNA_ID, lab), by="sgRNA_ID")

# ==============================================================================
# Plot growth curves, accounting for uncertainty
# ==============================================================================

# fig_growthcurve_R1_GSK_uncertainty <- plot_growthcurves_uncertainty(props_tidy_onoff, rdiff_df_uncertainty_GSK_summary, whichrep=1, whichdrug="GSK", sgRNA_ID_special="EZH2_030", bigval=bigpropcutoff, bigspots=c("switch"), whichweekmap=weekmap, whichguideattributes=guide_attributes)
# fig_growthcurve_R2_GSK_uncertainty <- plot_growthcurves_uncertainty(props_tidy_onoff, rdiff_df_uncertainty_GSK_summary, whichrep=2, whichdrug="GSK", sgRNA_ID_special="EZH2_030", bigval=bigpropcutoff, bigspots=c("switch"), whichweekmap=weekmap, whichguideattributes=guide_attributes)
# fig_growthcurve_R3_GSK_uncertainty <- plot_growthcurves_uncertainty(props_tidy_onoff, rdiff_df_uncertainty_GSK_summary, whichrep=3, whichdrug="GSK", sgRNA_ID_special="EZH2_030", bigval=bigpropcutoff, bigspots=c("switch"), whichweekmap=weekmap, whichguideattributes=guide_attributes)

# ggsave(fig_growthcurve_R1_GSK_uncertainty,file="figures/growthcurve_R1_GSK_uncertainty.pdf", width=2.5, height=2.5)
# ggsave(fig_growthcurve_R1_GSK_uncertainty,file="figures/growthcurve_R1_GSK_uncertainty.png", width=2.5, height=2.5)
# ggsave(fig_growthcurve_R2_GSK_uncertainty,file="figures/growthcurve_R2_GSK_uncertainty.pdf", width=2.5, height=2.5)
# ggsave(fig_growthcurve_R2_GSK_uncertainty,file="figures/growthcurve_R2_GSK_uncertainty.png", width=2.5, height=2.5)
# ggsave(fig_growthcurve_R3_GSK_uncertainty,file="figures/growthcurve_R3_GSK_uncertainty.pdf", width=2.5, height=2.5)
# ggsave(fig_growthcurve_R3_GSK_uncertainty,file="figures/growthcurve_R3_GSK_uncertainty.png", width=2.5, height=2.5)

# fig_growthcurve_R1_EED_uncertainty <- plot_growthcurves_uncertainty(props_tidy_onoff, rdiff_df_uncertainty_EED_summary, whichrep=1, whichdrug="EED", sgRNA_ID_special="EED_099", bigval=bigpropcutoff, bigspots=c("switch"), whichweekmap=weekmap, whichguideattributes=guide_attributes)
# fig_growthcurve_R2_EED_uncertainty <- plot_growthcurves_uncertainty(props_tidy_onoff, rdiff_df_uncertainty_EED_summary, whichrep=2, whichdrug="EED", sgRNA_ID_special="EED_099", bigval=bigpropcutoff, bigspots=c("switch"), whichweekmap=weekmap, whichguideattributes=guide_attributes)
# fig_growthcurve_R3_EED_uncertainty <- plot_growthcurves_uncertainty(props_tidy_onoff, rdiff_df_uncertainty_EED_summary, whichrep=3, whichdrug="EED", sgRNA_ID_special="EED_099", bigval=bigpropcutoff, bigspots=c("switch"), whichweekmap=weekmap, whichguideattributes=guide_attributes)

# ggsave(fig_growthcurve_R1_EED_uncertainty,file="figures/growthcurve_R1_EED_uncertainty.pdf", width=2.5, height=2.5)
# ggsave(fig_growthcurve_R1_EED_uncertainty,file="figures/growthcurve_R1_EED_uncertainty.png", width=2.5, height=2.5)
# ggsave(fig_growthcurve_R2_EED_uncertainty,file="figures/growthcurve_R2_EED_uncertainty.pdf", width=2.5, height=2.5)
# ggsave(fig_growthcurve_R2_EED_uncertainty,file="figures/growthcurve_R2_EED_uncertainty.png", width=2.5, height=2.5)
# ggsave(fig_growthcurve_R3_EED_uncertainty,file="figures/growthcurve_R3_EED_uncertainty.pdf", width=2.5, height=2.5)
# ggsave(fig_growthcurve_R3_EED_uncertainty,file="figures/growthcurve_R3_EED_uncertainty.png", width=2.5, height=2.5)


# ==============================================================================
# Plot scatter of proportion vs. rdiff
# ==============================================================================

rdiffpropdf_GSK <- rdiff_df_GSK %>% 
	# Append proportions: 
	left_join((props_tidy_onoff %>% 
			filter(drug=="GSK" & event=="switch") %>% 
			select(sgRNA_ID, drug, rep, prop)),
			by=c("sgRNA_ID","drug","rep")
		) %>% 
	# Deal with missing proportions for rep=-1 (late rdiff): 
	group_by(sgRNA_ID,drug) %>% 
	mutate(prop_avg=case_when(rep==0~prop,TRUE~-Inf)) %>% 
	mutate(prop_avg=max(prop_avg)) %>% 
	mutate(prop=case_when(rep==-1~prop_avg,TRUE~prop)) %>% 
	select(-prop_avg) %>% 
	# Join guide names: 
	left_join(select(guide_attributes,sgRNA_ID,lab),by="sgRNA_ID") %>% 
	mutate(lab=case_when(
		rep>0~paste0(lab,"\nR",rep),
		TRUE~paste0(lab,"\nAvg")))

# Extract just the replicates: 
rdiffpropdf_GSK_replicate <- rdiffpropdf_GSK %>% filter(rep>0)

# Extract the averages with early rdiff averaging: 
rdiffpropdf_GSK_avg <- rdiffpropdf_GSK %>% filter(rep==0) 

# Extract the averages with late rdiff averaging: 
rdiffpropdf_GSK_avg_laterdiff <- rdiffpropdf_GSK %>% filter(rep==-1) 


fig_GSK_scatter_avg <- ggplot() + 
	geom_hline(aes(yintercept=0), size=0.2, lty="dashed") + 
	geom_vline(aes(xintercept=bigpropcutoff), size=0.2, lty="dashed") + 
	geom_point(data=rdiffpropdf_GSK_avg, aes(x=prop,y=-rdiff), size=0.5, col="lightgray") + 
	geom_point(data=filter(rdiffpropdf_GSK_avg,rdiff<0 & prop>bigpropcutoff), aes(x=prop,y=-rdiff), size=1, col="red") + 
	geom_point(data=filter(rdiffpropdf_GSK_avg,rdiff>=0 & sgRNA_ID!=special_guide_GSK & prop>bigpropcutoff), aes(x=prop,y=-rdiff), size=1, col="blue") + 
	geom_point(data=filter(rdiffpropdf_GSK_avg, sgRNA_ID==special_guide_GSK), aes(x=prop,y=-rdiff), size=1, col="black") + 
	geom_text_repel(data=filter(rdiffpropdf_GSK_avg,rdiff<0 & prop>bigpropcutoff), aes(x=prop,y=-rdiff, label=lab), size=2, col="red", min.segment.length=0.1, segment.size=0.3, box.padding=0.2) + 
	geom_text_repel(data=filter(rdiffpropdf_GSK_avg,rdiff>=0 & sgRNA_ID!=special_guide_GSK & prop>bigpropcutoff), aes(x=prop,y=-rdiff, label=lab), size=2, col="blue", min.segment.length=0.1, segment.size=0.3, box.padding=0.2, alpha=0.5) + 
	geom_text_repel(data=filter(rdiffpropdf_GSK_avg,sgRNA_ID==special_guide_GSK), aes(x=prop,y=-rdiff, label=lab), size=2, col="black", min.segment.length=0.1, segment.size=0.3, box.padding=0.2, alpha=0.5) + 
	theme_classic() + 
	scale_x_continuous(breaks=seq(from=0,to=1,by=0.05), limits=c(0,0.3), expand=c(0.000,0.0001)) + 
	scale_y_continuous(breaks=seq(from=-10,to=10,by=5), limits=c(-10,10), expand=c(0.000,0.0001)) + 
	labs(title="GSK",x="Proportion of pool upon drug removal", y="Addiction score")

ggsave(fig_GSK_scatter_avg,file="figures/GSK_scatter_avg.pdf", width=2.5, height=2.5)
ggsave(fig_GSK_scatter_avg,file="figures/GSK_scatter_avg.png", width=2.5, height=2.5)

ggsave(fig_GSK_scatter_avg,file="figures/GSK_scatter_avg_rect.pdf", width=4, height=2.5)
ggsave(fig_GSK_scatter_avg,file="figures/GSK_scatter_avg_rect.png", width=4, height=2.5)

fig_GSK_scatter_avg_laterdiff <- ggplot() + 
	geom_hline(aes(yintercept=0), size=0.2, lty="dashed") + 
	geom_vline(aes(xintercept=bigpropcutoff), size=0.2, lty="dashed") + 
	geom_point(data=rdiffpropdf_GSK_avg_laterdiff, aes(x=prop,y=-rdiff), size=0.5, col="lightgray") + 
	geom_point(data=filter(rdiffpropdf_GSK_avg_laterdiff,rdiff<0 & prop>bigpropcutoff), aes(x=prop,y=-rdiff), size=1, col="red") + 
	geom_point(data=filter(rdiffpropdf_GSK_avg_laterdiff,rdiff>=0 & sgRNA_ID!=special_guide_GSK & prop>bigpropcutoff), aes(x=prop,y=-rdiff), size=1, col="blue") + 
	geom_point(data=filter(rdiffpropdf_GSK_avg_laterdiff, sgRNA_ID==special_guide_GSK), aes(x=prop,y=-rdiff), size=1, col="black") + 
	geom_text_repel(data=filter(rdiffpropdf_GSK_avg_laterdiff,rdiff<0 & prop>bigpropcutoff), aes(x=prop,y=-rdiff, label=lab), size=2, col="red", min.segment.length=0.1, segment.size=0.3, box.padding=0.2) + 
	geom_text_repel(data=filter(rdiffpropdf_GSK_avg_laterdiff,rdiff>=0 & sgRNA_ID!=special_guide_GSK & prop>bigpropcutoff), aes(x=prop,y=-rdiff, label=lab), size=2, col="blue", min.segment.length=0.1, segment.size=0.3, box.padding=0.2, alpha=0.5) + 
	geom_text_repel(data=filter(rdiffpropdf_GSK_avg_laterdiff,sgRNA_ID==special_guide_GSK), aes(x=prop,y=-rdiff, label=lab), size=2, col="black", min.segment.length=0.1, segment.size=0.3, box.padding=0.2, alpha=0.5) + 
	theme_classic() + 
	scale_x_continuous(breaks=seq(from=0,to=1,by=0.05)) + 
	labs(title="GSK",x="Proportion of pool upon drug removal", y="Addiction score")

ggsave(fig_GSK_scatter_avg_laterdiff,file="figures/GSK_scatter_avg_laterdiff.pdf", width=2.5, height=2.5)
ggsave(fig_GSK_scatter_avg_laterdiff,file="figures/GSK_scatter_avg_laterdiff.png", width=2.5, height=2.5)


fig_GSK_scatter_replicate <- ggplot() + 
	geom_hline(aes(yintercept=0), size=0.2, lty="dashed") + 
	geom_vline(aes(xintercept=bigpropcutoff), size=0.2, lty="dashed") + 
	geom_point(data=rdiffpropdf_GSK_replicate, aes(x=prop,y=-rdiff), size=0.5, col="lightgray") + 
	geom_point(data=filter(rdiffpropdf_GSK_replicate,rdiff<0 & prop>bigpropcutoff), aes(x=prop,y=-rdiff), size=1, col="red") + 
	geom_point(data=filter(rdiffpropdf_GSK_replicate,rdiff>=0 & sgRNA_ID!=special_guide_GSK & prop>bigpropcutoff), aes(x=prop,y=-rdiff), size=1, col="blue") + 
	geom_point(data=filter(rdiffpropdf_GSK_replicate, sgRNA_ID==special_guide_GSK), aes(x=prop,y=-rdiff), size=1, col="black") + 
	geom_text_repel(data=filter(rdiffpropdf_GSK_replicate,rdiff<0 & prop>bigpropcutoff), aes(x=prop,y=-rdiff, label=lab), size=2, col="red", min.segment.length=0.1, segment.size=0.3, box.padding=0.2) + 
	geom_text_repel(data=filter(rdiffpropdf_GSK_replicate,rdiff>=0 & sgRNA_ID!=special_guide_GSK & prop>bigpropcutoff), aes(x=prop,y=-rdiff, label=lab), size=2, col="blue", min.segment.length=0.1, segment.size=0.3, box.padding=0.2, alpha=0.5) + 
	geom_text_repel(data=filter(rdiffpropdf_GSK_replicate,sgRNA_ID==special_guide_GSK), aes(x=prop,y=-rdiff, label=lab), size=2, col="black", min.segment.length=0.1, segment.size=0.3, box.padding=0.2, alpha=0.5) + 
	theme_classic() + 
	scale_x_continuous(breaks=seq(from=0,to=1,by=0.05)) + 
	labs(title="GSK",x="Proportion of pool upon drug removal", y="Addiction score")

ggsave(fig_GSK_scatter_replicate,file="figures/GSK_scatter_replicate.pdf", width=2.5, height=2.5)
ggsave(fig_GSK_scatter_replicate,file="figures/GSK_scatter_replicate.png", width=2.5, height=2.5)


rdiffpropdf_EED <- rdiff_df_EED %>% 
	left_join((props_tidy_onoff %>% 
			filter(drug=="EED" & event=="switch") %>% 
			select(sgRNA_ID, drug, rep, prop)),
			by=c("sgRNA_ID","drug","rep")
		) %>% 
	# Deal with missing proportions for rep=-1 (late rdiff): 
	group_by(sgRNA_ID,drug) %>% 
	mutate(prop_avg=case_when(rep==0~prop,TRUE~-Inf)) %>% 
	mutate(prop_avg=max(prop_avg)) %>% 
	mutate(prop=case_when(rep==-1~prop_avg,TRUE~prop)) %>% 
	select(-prop_avg) %>% 
	left_join(select(guide_attributes,sgRNA_ID,lab),by="sgRNA_ID") %>% 
	mutate(lab=case_when(
		rep>0~paste0(lab,"\nR",rep),
		TRUE~paste0(lab,"\nAvg")))

rdiffpropdf_EED_replicate <- rdiffpropdf_EED %>% filter(rep>0)

rdiffpropdf_EED_avg <- rdiffpropdf_EED %>% filter(rep==0)

rdiffpropdf_EED_avg_laterdiff <- rdiffpropdf_EED %>% filter(rep==-1)

fig_EED_scatter_avg <- ggplot() + 
	geom_hline(aes(yintercept=0), size=0.2, lty="dashed") + 
	geom_vline(aes(xintercept=bigpropcutoff), size=0.2, lty="dashed") + 
	geom_point(data=rdiffpropdf_EED_avg, aes(x=prop,y=-rdiff), size=0.5, col="lightgray") + 
	geom_point(data=filter(rdiffpropdf_EED_avg,rdiff<0 & prop>bigpropcutoff), aes(x=prop,y=-rdiff), size=1, col="red") + 
	geom_point(data=filter(rdiffpropdf_EED_avg,rdiff>=0 & sgRNA_ID!=special_guide_EED & prop>bigpropcutoff), aes(x=prop,y=-rdiff), size=1, col="blue") + 
	geom_point(data=filter(rdiffpropdf_EED_avg, sgRNA_ID==special_guide_EED), aes(x=prop,y=-rdiff), size=1, col="black") + 
	geom_text_repel(data=filter(rdiffpropdf_EED_avg,rdiff<0 & prop>bigpropcutoff), aes(x=prop,y=-rdiff, label=lab), size=2, col="red", min.segment.length=0.1, segment.size=0.3, box.padding=0.2) + 
	geom_text_repel(data=filter(rdiffpropdf_EED_avg,rdiff>=0 & sgRNA_ID!=special_guide_EED & prop>bigpropcutoff), aes(x=prop,y=-rdiff, label=lab), size=2, col="blue", min.segment.length=0.1, segment.size=0.3, box.padding=0.2, alpha=0.5) + 
	geom_text_repel(data=filter(rdiffpropdf_EED_avg,sgRNA_ID==special_guide_EED), aes(x=prop,y=-rdiff, label=lab), size=2, col="black", min.segment.length=0.1, segment.size=0.3, box.padding=0.2, alpha=0.5) + 
	theme_classic() + 
	scale_x_continuous(breaks=seq(from=0,to=1,by=0.05), limits=c(0,0.3), expand=c(0.000,0.0001)) + 
	scale_y_continuous(breaks=seq(from=-10,to=10,by=5), limits=c(-10,10), expand=c(0.000,0.0001)) + 
	labs(title="EED",x="Proportion of pool upon drug removal", y="Addiction score")

ggsave(fig_EED_scatter_avg,file="figures/EED_scatter_avg.pdf", width=2.5, height=2.5)
ggsave(fig_EED_scatter_avg,file="figures/EED_scatter_avg.png", width=2.5, height=2.5)

ggsave(fig_EED_scatter_avg,file="figures/EED_scatter_avg_rect.pdf", width=4, height=2.5)
ggsave(fig_EED_scatter_avg,file="figures/EED_scatter_avg_rect.png", width=4, height=2.5)


fig_EED_scatter_avg_laterdiff <- ggplot() + 
	geom_hline(aes(yintercept=0), size=0.2, lty="dashed") + 
	geom_vline(aes(xintercept=bigpropcutoff), size=0.2, lty="dashed") + 
	geom_point(data=rdiffpropdf_EED_avg_laterdiff, aes(x=prop,y=-rdiff), size=0.5, col="lightgray") + 
	geom_point(data=filter(rdiffpropdf_EED_avg_laterdiff,rdiff<0 & prop>bigpropcutoff), aes(x=prop,y=-rdiff), size=1, col="red") + 
	geom_point(data=filter(rdiffpropdf_EED_avg_laterdiff,rdiff>=0 & sgRNA_ID!=special_guide_EED & prop>bigpropcutoff), aes(x=prop,y=-rdiff), size=1, col="blue") + 
	geom_point(data=filter(rdiffpropdf_EED_avg_laterdiff, sgRNA_ID==special_guide_EED), aes(x=prop,y=-rdiff), size=1, col="black") + 
	geom_text_repel(data=filter(rdiffpropdf_EED_avg_laterdiff,rdiff<0 & prop>bigpropcutoff), aes(x=prop,y=-rdiff, label=lab), size=2, col="red", min.segment.length=0.1, segment.size=0.3, box.padding=0.2) + 
	geom_text_repel(data=filter(rdiffpropdf_EED_avg_laterdiff,rdiff>=0 & sgRNA_ID!=special_guide_EED & prop>bigpropcutoff), aes(x=prop,y=-rdiff, label=lab), size=2, col="blue", min.segment.length=0.1, segment.size=0.3, box.padding=0.2, alpha=0.5) + 
	geom_text_repel(data=filter(rdiffpropdf_EED_avg_laterdiff,sgRNA_ID==special_guide_EED), aes(x=prop,y=-rdiff, label=lab), size=2, col="black", min.segment.length=0.1, segment.size=0.3, box.padding=0.2, alpha=0.5) + 
	theme_classic() + 
	scale_x_continuous(breaks=seq(from=0,to=1,by=0.05)) + 
	labs(title="EED",x="Proportion of pool upon drug removal", y="Addiction score")

ggsave(fig_EED_scatter_avg_laterdiff,file="figures/EED_scatter_avg_laterdiff.pdf", width=2.5, height=2.5)
ggsave(fig_EED_scatter_avg_laterdiff,file="figures/EED_scatter_avg_laterdiff.png", width=2.5, height=2.5)


fig_EED_scatter_replicate <- ggplot() + 
	geom_hline(aes(yintercept=0), size=0.2, lty="dashed") + 
	geom_vline(aes(xintercept=bigpropcutoff), size=0.2, lty="dashed") + 
	geom_point(data=rdiffpropdf_EED_replicate, aes(x=prop,y=-rdiff), size=0.5, col="lightgray") + 
	geom_point(data=filter(rdiffpropdf_EED_replicate,rdiff<0 & prop>bigpropcutoff), aes(x=prop,y=-rdiff), size=1, col="red") + 
	geom_point(data=filter(rdiffpropdf_EED_replicate,rdiff>=0 & sgRNA_ID!=special_guide_EED & prop>bigpropcutoff), aes(x=prop,y=-rdiff), size=1, col="blue") + 
	geom_point(data=filter(rdiffpropdf_EED_replicate, sgRNA_ID==special_guide_EED), aes(x=prop,y=-rdiff), size=1, col="black") + 
	geom_text_repel(data=filter(rdiffpropdf_EED_replicate,rdiff<0 & prop>bigpropcutoff), aes(x=prop,y=-rdiff, label=lab), size=2, col="red", min.segment.length=0.1, segment.size=0.3, box.padding=0.2) + 
	geom_text_repel(data=filter(rdiffpropdf_EED_replicate,rdiff>=0 & sgRNA_ID!=special_guide_EED & prop>bigpropcutoff), aes(x=prop,y=-rdiff, label=lab), size=2, col="blue", min.segment.length=0.1, segment.size=0.3, box.padding=0.2, alpha=0.5) + 
	geom_text_repel(data=filter(rdiffpropdf_EED_replicate,sgRNA_ID==special_guide_EED), aes(x=prop,y=-rdiff, label=lab), size=2, col="black", min.segment.length=0.1, segment.size=0.3, box.padding=0.2, alpha=0.5) + 
	theme_classic() + 
	scale_x_continuous(breaks=seq(from=0,to=1,by=0.05)) + 
	labs(title="EED",x="Proportion of pool upon drug removal", y="Addiction score")

ggsave(fig_EED_scatter_replicate,file="figures/EED_scatter_replicate.pdf", width=2.5, height=2.5)
ggsave(fig_EED_scatter_replicate,file="figures/EED_scatter_replicate.png", width=2.5, height=2.5)

# ==============================================================================
# Save addiction tables at various cutoffs
# ==============================================================================

rdifftab_output <- rbind(rdiffpropdf_GSK, rdiffpropdf_EED) %>% 
	select(-lab) %>% 
	left_join(select(guide_attributes, sgRNA_ID, lab), by="sgRNA_ID") 

rdifftab_wide_output_1 <- rdifftab_output %>% 
	mutate(addscore=1-rdiff) %>% 
	select(sgRNA_ID, drug, rep, addscore) %>% 
	mutate(rep=paste0("R",rep,"_addscore")) %>% 
	pivot_wider(names_from="rep",values_from="addscore")
 
rdifftab_wide_output_2 <- rdifftab_output %>% 
	select(sgRNA_ID, drug, rep, prop) %>% 
	mutate(rep=paste0("R",rep,"_prop")) %>% 
	pivot_wider(names_from="rep",values_from="prop")

rdifftab_wide_output <- inner_join(rdifftab_wide_output_1, rdifftab_wide_output_2, by=c("sgRNA_ID","drug")) %>% 
	select(sgRNA_ID, drug, 
		`R-1_addscore`, `R-1_prop`, 
		R0_addscore, R0_prop, 
		R1_addscore, R1_prop, 
		R2_addscore, R2_prop, 
		R3_addscore, R3_prop) %>% 
	left_join(select(guide_attributes, sgRNA_ID, Gene, sgRNA_seq, cut_site_AA, Domain), by="sgRNA_ID")

addictiontab_output <- rbind(rdiffpropdf_GSK, rdiffpropdf_EED) %>% 
	filter(rep>=0) %>% 
	filter(rdiff<0) %>% 
	filter(prop>0.01) %>% 
	mutate(propcat=case_when(prop>0.05~">5%",prop>0.025~">2.5%",TRUE~">1%")) %>%
	mutate(propcat=factor(propcat, levels=c(">1%",">2.5%",">5%"))) %>% 
	arrange(drug, desc(propcat), rdiff, rep) %>%
	select(-lab) %>% 
	left_join(select(guide_attributes, sgRNA_ID, lab), by="sgRNA_ID") 

write_csv(rdifftab_output, file="output/rdifftab.csv")
write_csv(rdifftab_wide_output, file="output/rdifftab_wide.csv")
write_csv(addictiontab_output, file="output/addictiontab.csv")

# ==============================================================================
# Plot histograms of rdiff by replicate
# ==============================================================================

fig_rdiffhist_GSK <- rdiff_df_GSK %>% 
	filter(rep>0) %>% 
	ggplot() + 
		annotate("rect", xmin=0, xmax=Inf, ymin=-Inf, ymax=Inf, fill="red", alpha=0.2) +
		geom_histogram(aes(x=-rdiff), binwidth=0.5, fill="lightblue", col="black", size=0.2) + 
		geom_vline(aes(xintercept=0), size=0.2) + 
		facet_wrap(~factor(rep)) + 
		theme_classic() + 
		labs(title="GSK", x="Addiction score", y="Number of guides")
ggsave(fig_rdiffhist_GSK,file="figures/rdiffhist_GSK.pdf", width=2.5, height=2.5)
ggsave(fig_rdiffhist_GSK,file="figures/rdiffhist_GSK.png", width=2.5, height=2.5)

fig_rdiffhist_EED <- rdiff_df_EED %>% 
	filter(rep>0) %>% 
	ggplot() + 
		annotate("rect", xmin=0, xmax=Inf, ymin=-Inf, ymax=Inf, fill="red", alpha=0.2) +
		geom_histogram(aes(x=-rdiff), binwidth=0.5, fill="lightblue", col="black", size=0.2) + 
		geom_vline(aes(xintercept=0), size=0.2) + 
		facet_wrap(~factor(rep)) + 
		theme_classic() + 
		labs(title="EED", x="Addiction score", y="Number of guides")
ggsave(fig_rdiffhist_EED,file="figures/rdiffhist_EED.pdf", width=2.5, height=2.5)
ggsave(fig_rdiffhist_EED,file="figures/rdiffhist_EED.png", width=2.5, height=2.5)

