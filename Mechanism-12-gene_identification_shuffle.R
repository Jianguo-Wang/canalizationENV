

### prediction lethality 
set.seed(123)
dataTemp_shuffle<-data.frame(
	x=M_D_GD_mae_x_Zscore,
	y=M_D_GD_mae_y_Zscore,
	z=M_D_GD_mae_z_Zscore,
	lethal=sample(embryo_lethal)
)


# Load required libraries
library(pROC)      # For ROC/AUC calculations
library(ggplot2)   # For enhanced plotting
library(caret)     # For data splitting (optional but recommended)
library(tidyverse) # For data manipulation

# ==================== STEP 1: Prepare Data ====================
# Use mtcars dataset, convert am to factor for binary classification

# Split data into training and testing sets (80/20 split)
set.seed(123)  # For reproducibility
train_index <- sample(1:dim(dataTemp_shuffle)[1],round(dim(dataTemp_shuffle)[1]*0.7))
train_data <- dataTemp_shuffle[train_index, ]
test_data <- dataTemp_shuffle[-train_index, ]

# ==================== STEP 2: Fit Logistic Regression Model ====================
# Model predicting transmission type (am) based on horsepower (hp) and weight (wt)
logit_model <- glm(as.numeric(lethal) ~ x+y+z, 
                   data = train_data, 
                   family = binomial(link = "logit"))

# View model summary
summary(logit_model)

# Interpretation of coefficients:
# - Coefficients are in log-odds
# - Positive coefficient = increases probability of am=1 (manual)
# - Negative coefficient = decreases probability of am=1

# ==================== STEP 3: Make Predictions ====================
# Get predicted probabilities for the test set
test_predictions <- predict(logit_model, 
                            newdata = test_data, 
                            type = "response")

# Convert to factor with 0.5 cutoff (standard threshold)
predicted_class <- ifelse(test_predictions > 0.5, 1, 0)
predicted_class <- factor(predicted_class, levels = c(0, 1))

# ==================== STEP 4: Calculate and Plot ROC/AUC ====================
# Calculate ROC curve
roc_result <- roc(response = as.numeric(test_data$lethal), 
                  predictor = test_predictions)

# Print AUC
cat("AUC:", auc(roc_result), "\n")

# Basic ROC plot
plot(roc_result, 
     main = paste("ROC Curve (AUC =", round(auc(roc_result), 3), ")"),
     col = "blue", 
     lwd = 2,
     print.auc = TRUE,
     print.auc.x = 0.6, 
     print.auc.y = 0.4)

# Add diagonal reference line
abline(a = 0, b = 1, lty = 2, col = "gray")

# ==================== OPTIONAL: Enhanced ROC Plot with ggplot2 ====================
# Create dataframe for ROC curve
roc_df_shuffle <- data.frame(
  sensitivity = roc_result$sensitivities,
  specificity = roc_result$specificities,
  thresholds = roc_result$thresholds
)

# Calculate 1 - specificity
roc_df_shuffle$one_minus_specificity <- 1 - roc_df_shuffle$specificity

# Create enhanced ROC plot
ggplot(roc_df_shuffle, aes(x = one_minus_specificity, y = sensitivity)) +
	geom_line(color = "blue", size = 0.5) +
	geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
	geom_polygon(aes(x = one_minus_specificity, y = sensitivity), fill = "skyblue", alpha = 0.3) +
	scale_x_continuous(breaks=c(0,0.5,1))+
	scale_y_continuous(breaks=c(0,0.5,1))+
	labs(x = "False Positive Rate",	y = "True Positive Rate", title = "Lethality prediction (shuffling phenotype)")+ #subtitle = paste("AUC =", round(auc(roc_result), 3)), 
	annotate("text", x = 0.7, y = 0.3, 
		   label = paste("AUC =", round(auc(roc_result), 3)),
		   size = 2, color = "darkred", fontface = "bold") 
	#theme(plot.title = element_text(hjust = 0.5, face = "regular"),
	#	plot.subtitle = element_text(hjust = 0.5))

ggsave("lethal_prediction_by_phi_model_shuffle.pdf",width=unit(2.1,"in"),height=unit(2.2,"in"))

