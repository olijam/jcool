ITEM_REPORT = $(HOME)/bin/item_report
BANFF_PARSER = $(HOME)/bin/banff_parser
BANFF_DIR = /data2/banff_by_day_and_campaign

all: day camp report

DAYS      = $(sort $(notdir $(wildcard $(BANFF_DIR)/*)))
DAY_FILES = $(foreach day, $(DAYS), day$(day).csv.gz)

day: $(DAY_FILES)

day%.csv.gz: $(BANFF_DIR)/% $(BANFF_PARSER)
	xzcat $(BANFF_DIR)/$*/*[0-9].xz | $(BANFF_PARSER) | gzip > $@

CAMPAIGNS  = $(sort $(basename $(notdir $(wildcard $(BANFF_DIR)/*/*[0-9].xz))))
CAMP_FILES = $(foreach camp, $(CAMPAIGNS), camp$(camp).csv.gz)

camp: $(CAMP_FILES)

camp%.csv.gz: $(wildcard $(BANFF_DIR)/*/%.xz) $(BANFF_PARSER)
	xzcat $(BANFF_DIR)/*/$*.xz | $(BANFF_PARSER) | gzip > $@
#	ls $(BANFF_DIR)/*/$*.xz
#

DAY_REPORT_FILES = $(foreach day, $(DAYS), day$(day)_report.csv)

day_report: $(DAY_REPORT_FILES)

CAMP_REPORT_FILES = $(foreach camp, $(CAMPAIGNS), camp$(camp)_report.csv)

camp_report: $(CAMP_REPORT_FILES) 

REPORT_FILES = $(DAY_REPORT_FILES) $(CAMP_REPORT_FILES)

report: $(REPORT_FILES)

%_report.csv : %.csv.gz $(ITEM_REPORT)
	$(ITEM_REPORT) $< > $@.log 2>&1

clean:
	find ./ -size 20c -exec {} \;
