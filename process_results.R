suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("purrr"))
suppressPackageStartupMessages(library("ggplot2"))

my_path <-
  if (interactive()) {
    "."
  } else {
    if (length(commandArgs(TRUE)))
      commandArgs(TRUE)[1]
    else
      "."
  }

getdata <- function(x) {
  readRDS(file.path(my_path, x)) %>%
  select(G.r.p_A, G.r.p_B, G.r.p_AB) %>%
  pivot_longer(G.r.p_A:G.r.p_AB, names_to = "version") %>%
  separate(version, c("model", "rb", "effect"), "\\.")
}

allfiles <- dir(my_path, "^ckm.*\\.rds$")

if (length(allfiles) == 0L) {
  stop("No output files found in ",
       if (my_path == ".") "current path" else my_path)
}

message("Processing ", length(allfiles), " output files...")

alldat <- tibble(fname = allfiles) %>%
  separate(fname, c("j", "nmc", "ns", "nobs", "A", "B", "AB",
		    "ri0", "ri1", "rs0", "rs1", "version",
		    "date", "host", "pid"), "_",
	   remove = FALSE, convert = TRUE) %>%
  mutate(dat = map(fname, getdata)) %>%
  unnest(c(dat)) %>%
  mutate(sig = value < .05,
	 factor = factor(sub("^p_", "", effect), levels = c("A", "B", "AB")),
	 model = factor(if_else(model == "G", "GAMM", "LMEM")),
	 effsize = as.numeric(case_when(factor == "A" ~ A,
					factor == "B" ~ B,
					TRUE ~ AB)),
	 allocation = if_else(rb == "r", "randomized", "blocked"))

cases <- tibble(version = 1:8,
		case = factor(c("1. varying phase",
				"2. varying amp",
				"3. random walk 1",
				"4. random walk 2",
				"5. multi 1 + 3",
				"6. multi 1 + 4",
				"7. multi 2 + 3",
				"8. multi 2 + 4")))

mstats <- alldat %>%
  filter(version != 0L) %>%
  inner_join(cases, "version") %>%
  group_by(ns, nobs, effsize, case, model, allocation, factor) %>%
  summarize(psig = mean(sig), N = n()) %>%
  ungroup() %>%
  arrange(allocation, effsize, case, factor, model) %>%
  select(-ns, -nobs) %>%
  mutate(model = "GAMM-regfs")

sres <- readRDS("../data_derived/simulation_results.rds")

orig <- sres %>%
  filter(case %in% c(1L, 4L)) %>%
  select(case, A, B, AB,
         G.r.p_A, G.r.p_B, G.r.p_AB,
         L.r.p_A, L.r.p_B, L.r.p_AB) %>%
  pivot_longer(cols = c(G.r.p_A, G.r.p_B, G.r.p_AB,
                        L.r.p_A, L.r.p_B, L.r.p_AB),
               names_to = "version") %>%
  separate(version, c("model", "rb", "effect"), "\\.") %>%
  mutate(sig = value < .05,
	 factor = factor(sub("^p_", "", effect), levels = c("A", "B", "AB")),
	 model = factor(if_else(model == "G", "GAMM", "LMEM")),
	 effsize = as.numeric(case_when(factor == "A" ~ A,
					factor == "B" ~ B,
					TRUE ~ AB)),
	 allocation = if_else(rb == "r", "randomized", "blocked"))

orig_stats <- orig %>%
  rename(version = case) %>%
  inner_join(cases, "version") %>%
  group_by(effsize, case, model, allocation, factor) %>%
  summarize(psig = mean(sig), N = n()) %>%
  ungroup() %>%
  arrange(allocation, effsize, case, factor, model) 

result <- bind_rows(orig_stats, mstats)

effA <- ggplot(
  result %>% filter(factor == "A"),
  aes(effsize, psig)) +
  geom_point(aes(color = model, shape = model), alpha = .5) +
  geom_line(aes(color = model), alpha = .5) +
  scale_x_continuous(breaks = seq(0, .25, .05)) +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, .25)) +
  facet_grid(allocation ~ case) +
  theme(legend.position = "top",
	legend.margin = margin(0, 0, 0, 0),
	legend.box.margin = margin(-5, 0, -5, 0),        
	axis.text.x = element_text(size = 5),
	axis.text.y = element_text(size = 7)) +
  labs(x = "raw effect", y = "proportion significant")

effB <- ggplot(
  result %>% filter(factor == "B"),
  aes(effsize, psig)) +
  geom_point(aes(color = model, shape = model), alpha = .5) +
  geom_line(aes(color = model), alpha = .5) +
  scale_x_continuous(breaks = seq(0, .5, .1)) +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, .5)) +  
  facet_grid(allocation ~ case) +
  theme(legend.position = "top",
	legend.margin = margin(0, 0, 0, 0),
	legend.box.margin = margin(-5, 0, -5, 0),        
	axis.text.x = element_text(size = 7),
	axis.text.y = element_text(size = 7)) +
  labs(x = "raw effect", y = "proportion significant")

effAB <- ggplot(
  result %>% filter(factor == "AB"),
  aes(effsize, psig)) +
  geom_point(aes(color = model, shape = model), alpha = .5) +
  geom_line(aes(color = model), alpha = .5) +
  scale_x_continuous(breaks = seq(0, .5, .1)) +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, .5)) +  
  facet_grid(allocation ~ case) +
  theme(legend.position = "top",
	legend.margin = margin(0, 0, 0, 0),
	legend.box.margin = margin(-5, 0, -5, 0),        
	axis.text.x = element_text(size = 7),
	axis.text.y = element_text(size = 7)) +
  labs(x = "raw effect", y = "proportion significant")

ggsave("eff_A.png", effA, width = 10, height = 5)
ggsave("eff_B.png", effB, width = 10, height = 5)
ggsave("eff_AB.png", effAB, width = 10, height = 5)
message("wrote files 'eff_A.png', 'eff_B.png', 'eff_AB.png'")
