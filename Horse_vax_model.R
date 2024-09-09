################################################################################
#                                                                              #
# This code simulates the population of horses in Canada and by province, and  # 
# estimates the proportion of horses vaccinated for West Nile virus.           #
#                                                                              #
# Author: Pascale Aubry (with help from Josh Persi)                            #
#                                                                              #
################################################################################

#Load required packages
library(here)
library(mc2d)
library(tidyverse)

#Load the distribution parameters by province
dist_pars <- 
  read.csv(here::here("Horse_vax_data.csv"))

#The proportion of horses vaccinated for the first time is represented by a pert distribution
# The same distribution is repeated for all provinces.
primo_min <- 0
primo_max <- 0.1

#######################################################################################################

#The distribution of the number of doses sold/administered, and the distribution
# of the horse population are generated for each province.

# This is done using purrr:map() to iterate over the number of iterations, and for each
# iteration, pulling the appropriate value from the appropriate column. This is
# done in a rowwise fashion thanks to the .by argument to mutate, where we are
# in essence saying run the map function for each region.

#The number of horses vaccinated is estimated from the number of doses sold/administered,
# and the probability of primo-vaccination, knowing that 2 doses are administered in primo
# vaccination and 1 dose for booster vaccination
#A binomial distribution is used with n = number of doses sold/administered,
# and p = 1 / (2*primo + 1*(1-primo))

set.seed(1234)

n_iter <- 100

#The probability of primo vaccination is the same for all provinces
prob_primo <- 
  runif(n_iter, primo_min, primo_max)

#The number of doses and horse population is estimated for each province
dist_pars <- 
  dist_pars %>%
  mutate(
    output = list(map(
      .x = 1:n_iter,
      .f = function(x) {
        
        doses <- round(
          rpert(1, doses_min, doses_ML, doses_max),
          0
        )
        
        horse_pop <- round(
          rpert(1, pop_min, pop_ML, pop_max),
          0
        )
        
        output <- cbind.data.frame(
          iter = x,
          region = region,
          doses = doses,
          horse_pop = horse_pop
        )
      }
    )),
    .by = region
  )

# We now have a list column with each region containing a list equal to the
# number of iterations. Merge the list elements, which are dataframes, together.
dist_pars <- 
  dist_pars %>% 
  mutate(
    output = map(output, list_c),
    .by = region
  )

# We can unnest the nested dataframes of output 
output <- 
  dist_pars  %>% 
  select(output)  %>% 
  unnest(output)

#Now is time to add the probability primo vaccination (same values for all provinces)
output <- 
  output %>% 
  mutate(
    primo = rep(prob_primo, nrow(dist_pars))) %>% 
  #and simulate number horses vaccinated
  rowwise() %>% 
  mutate(
    n_vax = rbinom(1, doses, 1 / (2*primo + 1*(1-primo))),
    prop_vax = n_vax / horse_pop
  )

#Add up the provincial horse populations and number of vaccinated horses,
#by iteration, to obtain country-level estimates
output_Canada <- 
  output %>% 
  group_by(iter) %>% 
  summarise(
    horse_pop=sum(horse_pop),
    n_vax =sum(n_vax )) %>% 
  mutate(region = "Canada",
         prop_vax = n_vax / horse_pop) %>% 
  select(iter, region,  horse_pop, n_vax, prop_vax)

results <-
  bind_rows(output,
            output_Canada) %>% 
  select(-doses, -primo)

# Summary function to obtain median, mean, 95% and 99% confidence interval
my_summary <- function(arg1)  {
  min <- min(arg1, na.rm=TRUE)
  P0.5 <- quantile(arg1, probs=.005, na.rm=TRUE)
  P2.5 <- quantile(arg1, probs=.0255, na.rm=TRUE)
  med <- median(arg1, na.rm=TRUE)
  m <- mean(arg1, na.rm=TRUE)
  P97.5 <- quantile(arg1, probs=.975, na.rm=TRUE)
  P99.5 <- quantile(arg1, probs=.995, na.rm=TRUE)
  max <- max(arg1, na.rm=TRUE)
  data.frame(Min.=min, P0.5=P0.5, P2.5 =P2.5 , Med.=med, 
             Mean = m, P97.5=P97.5, P99.5=P99.5, Max.=max) 
}

#Summary stats for proportion vaccinated, by region
results %>% 
  group_by(region)   %>% 
  summarise(as_tibble(rbind(my_summary(prop_vax))))

#Summary stats for horse population, by region
results %>% 
  group_by(region)   %>% 
  summarise(as_tibble(rbind(my_summary(horse_pop))))


#Get the number of unique pairwise combinations of regions
n_pairs <-
  length(
    combn(
      dist_pars$region, 
      m=2, 
      simplify = FALSE)
  )

# Obtain all the pairwise differences in proportion vaccinated (between provinces)
prop_vax_diff <- results %>%
  filter(region != "Canada") %>% 
  # Nest and fully cross the data
  nest(data = prop_vax, .by = region) %>%
  expand(
    nesting(reference_region = region, reference_data = data),
    nesting(comparison_region = region, comparison_data = data)
  ) %>%
  # Remove records where the reference region and the comparison region are the
  # same (i.e. diff == 0)
  filter(reference_region != comparison_region) %>%
  # From each value of prop_vax per reference region, subtract each value
  # of prop_vax per comparison region
  mutate(difference = map2(reference_data, comparison_data, \(x, y) x - y)) %>%
  # Unnest the difference data
  unnest(difference) %>%
  # Select the relevant columns and rename the prop_vax column to show
  # it is a difference
  select(reference_region, comparison_region, prop_vax_diff = prop_vax) %>%
  # Group data by row
  rowwise() %>%
  # Per row, combine reference and comparison regions, sort them alphabetically,
  # and concatenate them so AB-Atlantic and Atlantic-AB comparisons both read
  # AB-Atlantic
  mutate(
    comparison_id = 
      list(c(reference_region, comparison_region)) %>%
      map_chr(
        \(x) x %>%
          str_sort() %>%
          str_c(collapse = "-")
      )
  ) %>%
  # Ungroup the dataset
  ungroup() %>%
  # Add the iteration number
  mutate(iter = rep(1:n_iter, times = 2 *n_pairs)) %>% 
  # Only keep unique comparison IDs
  group_by(iter,comparison_id) %>% 
  distinct(comparison_id, .keep_all = TRUE)

#Summary stats for pairwise differences in proportion vaccinated
summary_prop_vax_diff <-
  prop_vax_diff %>% 
  group_by(comparison_id)   %>% 
  summarise(
    as_tibble(
      rbind(
        my_summary(prop_vax_diff))))

#Find the pairs for which zero is NOT included in the interval [P0.5, P99.5] 
#These are the ones for which the prop_vax difference is statistically significant  
summary_prop_vax_diff %>% 
  filter((P0.5 * P99.5) > 0)
