data_first_inf_variant <- VASCO_data_org |> 
  mutate(Variant = Variant |> factor())
# 1. Build model
# model with custom knot locations
knots_loc <- c(seq(-300,-50, 30), -30,-20, seq(-10,30,5), 40,50, seq(70, 300, 30))

## fitting nc spline with fixed positions
cr_fit_fix_variant <- gam(CIS_subj_fat ~ s(Time, bs = 'cr', k = length(knots_loc), by = Variant) + Variant,
                  knots = list(Time = knots_loc),
                  data = data_first_inf_variant)
# CIS_subj_fat = CIS score
# Time = tijd sinds positieve test (numeric)
# Variant = variant van infectie (Delta of Omicron)

# 2. Predict
time_range <- range(data_first_inf_variant$Time)
variant_range <- levels(data_first_inf_variant$Variant)
data_frame_variant <- expand.grid(seq(from = time_range[1], to = time_range[2],1),
                                  variant_range)  
colnames(data_frame_variant) <- c("Time", "Variant")
data_frame_variant <- data.frame(data_frame_variant) 


pred_variant <- predict(cr_fit_fix_variant, newdata=data_frame_variant, se.fit=TRUE, discrete=F,
                      newdata.guaranteed=TRUE)
data_frame_variant <- data_frame_variant %>%
  mutate(prediction = pred_variant$fit,
         se = pred_variant$se.fit,
         ll = pred_variant$fit-1.96*pred_variant$se.fit,
         ul = pred_variant$fit+1.96*pred_variant$se.fit)

# 3. Plot
plot_cr_fit_fix_variant <- ggplot() +
  geom_line(data = data_frame_variant, aes(x=Time, y=prediction, colour = Variant)) +
  geom_ribbon(data = data_frame_variant,
              aes(x = Time, ymin = ll, ymax = ul, fill=Variant),
              alpha = 0.2)  +
  scale_x_continuous(breaks=seq(-300, 300, by=50)) +
  geom_vline(aes(xintercept = c(0,90,180,270)),linetype="dashed", linewidth=0.3, colour="red") +
  labs(title = "natural cubic spline fixed knots") +
  ylim(c(18, 30)) +
  theme_bw() +
  theme(
    axis.title = element_text(size=15),
    axis.text = element_text(size=14),
    strip.text = element_text(size=15),
    legend.text = element_text(size=14),
    legend.title = element_text(size=15),
  ); plot_cr_fit_fix_variant

# prediction at 3 months pre-infection, and 3, 6, 9 months post-infection
data_frame_variant[data_frame_variant$Time==-90,]
data_frame_variant[data_frame_variant$Time==90,]
data_frame_variant[data_frame_variant$Time==180,]
data_frame_variant[data_frame_variant$Time==270,]

# check op die tijdspunten of er significant verschil is in score tussen Delta en Omicron
diff_variant <- data_frame_variant %>% 
  filter(Time==-90 | Time==90 | Time==180 | Time == 270) %>% 
  select(Time, Variant, prediction, se, ll, ul) %>%
  arrange(Time) %>%
  group_by(Time) %>%
  mutate(difference = abs(prediction - lag(prediction))) %>%
  mutate(difference_se = sqrt((se^2)+(lag(se^2)))) %>%
  mutate(ll_perc = difference-1.96*difference_se,
         ul_perc = difference+1.96*difference_se) %>% ungroup()

