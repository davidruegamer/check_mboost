#### Iris example ####
data(iris)
### glmboost with multiple variables and intercept
iris$setosa <- factor(iris$Species == "setosa")
iris_glm <- glmboost(setosa ~ 1 + Sepal.Width + Sepal.Length + Petal.Width +
                       Petal.Length,
                     data = iris, control = boost_control(mstop = 50), 
                     family = Binomial(link = c("logit")))

# try extract_cum_expl_risk
expl_risk <- extract_cum_expl_risk(iris_glm)
plot(expl_risk, col = selected(iris_glm))

# try extract_selection_path
sel_path <- extract_sel_path(iris_glm)
barplot(as.numeric(sel_path), col = selected(iris_glm))