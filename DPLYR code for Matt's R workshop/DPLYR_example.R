library("dplyr", "tidyr")

setwd("O:/Coursework/DPLYR code for Matt's R workshop/")

# Genotype test: PCR sequence of partial HIV pol gene. Patient sequence compared to database to identify resistance causing mutations. 
# The 3 classes of drugs tested are NRTI, NNRTI, PI. 
# For this dataset, no drug resistance for that class = 0, resistance = 1

# Load database results for patients in study
db <- read.table("DB_results.csv", header=TRUE, sep=",", stringsAsFactors = FALSE, strip.white=TRUE)
attach(db)

# Load table with all the sequences for the clinic
seqs <- read.table("All_sequences.csv", header=TRUE, sep=",", stringsAsFactors = FALSE, strip.white=TRUE)
attach(seqs)


# We want to merge the database results with the sequences and create a few more variables to continue analysis
# Left join = only the data matching our db results will be added to the table 
  # (i.e. if the patient is not in our study, we don't want their sequence data)
  # Join by "ENUM" column because each genotype has unique ENUM identifier, there are multiple tests/patient so that is not unique

all <- left_join(db, seqs, by = "ENUM")

# look at table, notice there are mulitple columns with same data
head(all)

# get rid of duplicate columns and clean up column names
all_clean <- select(all, PATID.x, First_pos_test, ENUM, DRAWN_DATE.x, NRTI, NNRTI, PI, LOCATION, STRING)
colnames(all_clean) <- c("Pat_ID", "First_positive_test", "Geno_ID", "Geno_date", "NRTI", "NNRTI", "PI", "Location", "Sequence")

# Manipulate the table 
meta <-
  select(all_clean, -Sequence) %>%  # remove sequences for now        
  select(-NRTI, -NNRTI, -PI, everything()) %>%  # reorder remaining columns
  mutate(sum_resistance=NRTI+NNRTI+PI) # make new column of the total resistance for each genotype result
  
# Currently, HIV is not cured by the ARV drugs, only suppressed in the body. Therefore patients can cumulate resistance over time.
# If you look at the data, you can see that some patients have different resistance results. 
# This is because the PCR used in the genotype only records the most abundant HIV variant circulating at that time but there are
# many variants in the population of that patient so we often assume that resistance is cumulative and therefore we want to 
# categorize patients by their maximum resistance score ever (i.e. cumulative resistance).
  
cumulative_res <-
  meta %>%
  group_by(Pat_ID) %>%
  summarize(max_res = max(sum_resistance)) %>%
  select(Pat_ID, max_res) %>%
  arrange(max_res)  
  
# what is the number of patients in each resistant category?
count(cumulative_res, max_res)

# Add cumulative resistance category to the metadata
meta <- left_join(meta, cumulative_res, by="Pat_ID")
  
# Count the patients in each category
ccr <- select(meta, Pat_ID, max_res)
ccr_counts <- ccr %>% 
  group_by(Pat_ID,max_res) %>%
  summarise(freq = length(max_res)) 

# Basic statistical summary by cumulative resistance class
bymaxclass <- ccr_counts %>% group_by(max_res) %>% select(max_res, freq)
tapply(bymaxclass$freq, bymaxclass$max_res, summary)  
  
# Plot the results in a boxplot
boxplot(ccr_counts$freq ~ ccr_counts$max_res, xlab="Cumulative class resistance/patient", 
        ylab="frequency of tests", main="Number of genotypes/patient by cumulative class resistance")

