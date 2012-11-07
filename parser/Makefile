BANFF_DIR = /data2/banff_by_day_and_campaign

all: day camp report

DAYS      = $(sort $(notdir $(wildcard $(BANFF_DIR)/*)))
DAY_FILES = $(foreach day, $(DAYS), day$(day).csv.gz)

day: $(DAY_FILES)

day%.csv.gz: $(BANFF_DIR)/%
	xzcat $(BANFF_DIR)/$*/*[0-9].xz | banff_parser | gzip > $@

CAMPAIGNS  = $(sort $(basename $(notdir $(wildcard $(BANFF_DIR)/*/*[0-9].xz))))
CAMP_FILES = $(foreach camp, $(CAMPAIGNS), camp$(camp).csv.gz)

camp: $(CAMP_FILES)

camp%.csv.gz: $(wildcard $(BANFF_DIR)/*/%.xz)
	xzcat $(BANFF_DIR)/*/$*.xz | banff_parser | gzip > $@
#	ls $(BANFF_DIR)/*/$*.xz
#

DAY_REPORT_FILES = $(foreach day, $(DAYS), day$(day)_report.csv)

day_report: $(DAY_REPORT_FILES)

CAMP_REPORT_FILES = $(foreach camp, $(CAMPAIGNS), camp$(camp)_report.csv)

camp_report: $(CAMP_REPORT_FILES)

REPORT_FILES = $(DAY_REPORT_FILES) $(CAMP_REPORT_FILES)

report: $(REPORT_FILES)

%_report.csv : %.csv.gz
	item_report $< > $@.log 2>&1
