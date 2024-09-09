################################################################################
#                                                                              #
# This code performs a sensitivity analysis on the model to estimate the       # 
# proportion of horses vaccinated for West Nile virus in Canada.               #
#                                                                              #
# Author: Pascale Aubry                                                        #
#                                                                              #
################################################################################


### This code reproduces the left-hand plot of Figure 5.21 in Vose D, 2008. Risk analysis: a quantitative guide, 3rd ed. 735 pp. Chichester, England: Wiley.

# It is a tornado chart will all the input variables for the vaccination coverage in Canada
# It is the crudest type of sensitivity analysis, showing the Spearman rank 
# correlation coefficients between all the inputs and the vaccination coverage for Canada

# It gives an idea of the most important input variables (those that have at least a quarter of the maximum observed correlation)

#Also scatter plots of inputs vs ouput, and boxplot of vaccine coverage by province

#First run the vaccination model
source(here::here("FINAL horse vaccination model.R"))


# Correlations using all the provincial variables (populations and doses)-----------------
corr_data <- 
  output %>% 
  pivot_longer(
    cols = c("doses", "horse_pop"),
    names_to = "variable",
    values_to = "value"
  ) %>% 
  mutate(
    variable = paste(region, "_",variable, sep = "")  
  ) %>%
  select(iter, variable, value) 

corr_data <-
  output_Canada %>% 
  pivot_longer(
    cols = prop_vax,
    names_to = "variable",
    values_to = "value") %>% 
  select(iter, variable, value) %>% 
  bind_rows(corr_data, .) %>% 
  arrange(iter)

corr_data_wide <-
  corr_data %>% 
  pivot_wider(
    names_from = variable,
    values_from = value
  ) %>% 
  mutate(primo = prob_primo) %>% 
  select(prop_vax, everything())

all_corr <-
  round(cor(corr_data_wide,
            method="spearman"),
        digits = 3 # rounded to 2 decimals
  )

correlations <-
  enframe(all_corr[3:17,"prop_vax"], name = "variable", value = "coef") %>% 
  arrange(desc(abs(coef)))

#Order plot with stronger correlations on top, using the absolute value of the coefficient
tornado_plot_supp_mat <-
  ggplot(
    correlations, 
    aes(x = reorder(variable, abs(coef)), y = coef)) +
  geom_bar(stat="identity", color='blue',fill='blue') +
  coord_flip() + 
  labs(
    y = "Spearman's rank correlation coefficient",
    x = "")

tiff(file="tornado_plot2.tiff",
     width=6, height=4, units="in", res=300)
tornado_plot_supp_mat
dev.off()


### Have a look at what happens at the provincial level--------------------------

#Manitoba
corr_MB <- 
  output %>% 
  filter(region == "MB") %>% 
  pivot_longer(
    cols = c("doses", "horse_pop", "primo"),
    names_to = "variable",
    values_to = "value"
  )  %>%
  select(iter, variable, value) 

corr_MB <-
  output %>% 
  filter(region == "MB") %>% 
  pivot_longer(
    cols = prop_vax,
    names_to = "variable",
    values_to = "value") %>% 
  select(iter, variable, value) %>% 
  bind_rows(corr_MB, .) %>% 
  arrange(iter)

corr_MB_wide <-
  corr_MB %>% 
  pivot_wider(
    names_from = variable,
    values_from = value
  ) %>% 
  select(-iter)

MB_corr <-
  round(cor(corr_MB_wide,
            method="spearman"),
        digits = 3 # rounded to 2 decimals
  )

MB_correlations <-
  enframe(MB_corr[1:3,"prop_vax"], name = "variable", value = "coef") %>% 
  arrange(desc(abs(coef)))

#Order plot with stronger correlations on top, using the absolute value of the coefficient
ggplot(
  MB_correlations, 
  aes(x = reorder(variable, abs(coef)), y = coef)) +
  geom_bar(stat="identity", color='blue',fill='blue') +
  coord_flip()

#Ontario
corr_ON <- 
  output %>% 
  filter(region == "ON") %>% 
  pivot_longer(
    cols = c("doses", "horse_pop", "primo"),
    names_to = "variable",
    values_to = "value"
  )  %>%
  select(iter, variable, value) 

corr_ON <-
  output %>% 
  filter(region == "ON") %>% 
  pivot_longer(
    cols = prop_vax,
    names_to = "variable",
    values_to = "value") %>% 
  select(iter, variable, value) %>% 
  bind_rows(corr_ON, .) %>% 
  arrange(iter)

corr_ON_wide <-
  corr_ON %>% 
  pivot_wider(
    names_from = variable,
    values_from = value
  ) %>% 
  select(-iter)

ON_corr <-
  round(cor(corr_ON_wide,
            method="spearman"),
        digits = 3 # rounded to 2 decimals
  )

ON_correlations <-
  enframe(ON_corr[1:3,"prop_vax"], name = "variable", value = "coef") %>% 
  arrange(desc(abs(coef)))

#Order plot with stronger correlations on top, using the absolute value of the coefficient
ggplot(
  ON_correlations, 
  aes(x = reorder(variable, abs(coef)), y = coef)) +
  geom_bar(stat="identity", color='blue',fill='blue') +
  coord_flip()


# Scatter plots-----------------------------------------------------------------
inputs <-
  list_c(dist_pars$output)

sens <-
  right_join(
    inputs, results, by = c("iter", "region", "doses", "horse_pop")) 

sens_Can <-
  sens %>%
  filter(region == "Canada")

sens_provs <-
  sens %>%
  filter(region != "Canada")

#Scatter plots by province
sens_provs %>% 
  ggplot(aes(x=horse_pop, y=prop_vax)) + 
  geom_point(aes(color = factor(region)))

sens_provs %>% 
  ggplot( aes(x=doses, y=prop_vax)) + 
  geom_point(aes(color = factor(region)))

sens_provs %>% 
  ggplot( aes(x=primo, y=prop_vax)) + 
  geom_point(aes(color = factor(region)))

#Box plot vaccination coverage by province
sens_provs %>% 
  ggplot(aes(region, prop_vax)) +
  geom_boxplot()

