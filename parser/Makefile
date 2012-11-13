ZIPPER = /usr/bin/xz -4
ITEM_REPORT = $(HOME)/bin/item_report
BANFF_PARSER = $(HOME)/bin/banff_parser
BANFF_DIR = /data2/banff_by_day_and_campaign
LOG_DIR = log

.PHONY: all
all: day camp report

DAY_DIR   = days
DAYS      = $(sort $(notdir $(wildcard $(BANFF_DIR)/*)))
DAY_FILES = $(foreach day, $(DAYS), $(DAY_DIR)/day$(day).csv.xz)

.PHONY: day
day: $(DAY_FILES)
$(DAY_DIR)/day%.csv.xz: $(wildcard $(BANFF_DIR)/%/*[0-9].xz) $(BANFF_PARSER)
	xzcat $(BANFF_DIR)/$*/*[0-9].xz | $(BANFF_PARSER) | $(ZIPPER) > $@

CAMPAIGN_DIR = campaigns
CAMPAIGNS  = $(sort $(basename $(notdir $(wildcard $(BANFF_DIR)/*/*[0-9].xz))))
CAMP_FILES = $(foreach camp, $(CAMPAIGNS), $(CAMPAIGN_DIR)/camp$(camp).csv.xz)

.PHONY: camp
camp: $(CAMP_FILES)

vpath %.xz $(wildcard $(BANFF_DIR)/*)
$(CAMPAIGN_DIR)/camp%.csv.xz: %.xz $(BANFF_PARSER)
	xzcat $(BANFF_DIR)/*/$*.xz | $(BANFF_PARSER) | $(ZIPPER) > $@

ITEM_DIR = item_reports

.PHONY: day_report
DAY_REPORT_FILES = $(foreach day, $(DAYS), $(ITEM_DIR)/day$(day)_report.csv)
day_report: $(DAY_REPORT_FILES)
vpath day%.csv.xz $(DAY_DIR)
$(ITEM_DIR)/day%_report.csv: day%.csv.xz $(ITEM_REPORT)
	$(ITEM_REPORT) $< > $(LOG_DIR)/day$*_report.log 2>&1

.PHONY: camp_report
CAMP_REPORT_FILES = $(foreach camp, $(CAMPAIGNS), $(ITEM_DIR)/camp$(camp)_report.csv)
camp_report: $(CAMP_REPORT_FILES) 
vpath camp%.csv.xz $(CAMPAIGN_DIR)
$(ITEM_DIR)/camp%_report.csv: camp%.csv.xz $(ITEM_REPORT)
	$(ITEM_REPORT) $< > $(LOG_DIR)/camp$*_report.log 2>&1

.PHONY: report
REPORT_FILES = $(DAY_REPORT_FILES) $(CAMP_REPORT_FILES)
report: $(REPORT_FILES)

.PHONY: clean
clean:
	find ./ -name '*gz' -size 20c -ls -exec rm {} \;
	find ./ -name '*xz' -size 60c -ls -exec rm {} \;
	find ./ -name '*csv' -size 0c -ls -exec rm {} \;
