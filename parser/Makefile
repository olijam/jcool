BANFF_DIR = /data2/banff_by_day_and_campaign

CAMPAIGNS  = $(sort $(basename $(notdir $(wildcard $(BANFF_DIR)/*/*[0-9].xz))))
CAMP_FILES = $(foreach camp, $(CAMPAIGNS), camp$(camp).csv.gz)

DAYS      = $(sort $(notdir $(wildcard $(BANFF_DIR)/*)))
DAY_FILES = $(foreach day, $(DAYS), day$(day).csv.gz)

all: day camp

day: $(DAY_FILES)

day%.csv.gz: $(BANFF_DIR)/%
	xzcat $(BANFF_DIR)/$*/*[0-9].xz | banff_parser | gzip > $@

camp: $(CAMP_FILES)

camp%.csv.gz: $(wildcard $(BANFF_DIR)/*/%.xz)
	xzcat $(BANFF_DIR)/*/$*.xz | banff_parser | gzip > $@
#	ls $(BANFF_DIR)/*/$*.xz
#

report: *_report.c
