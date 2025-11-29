import json
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# --- Configuration ---
JSON_FILE = '..\..\model_results\protein_analysis_metrics.json'
FEATURE_IMPORTANCE_KEY = 'feature_importance'
TOP_N_FEATURES = 5 # You can change this to show more or fewer features

def load_and_process_data(file_path):
    """
    Loads the JSON data and restructures the feature importance metrics
    into a long-format Pandas DataFrame suitable for visualization.
    """
    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
    except FileNotFoundError:
        print(f"Error: JSON file not found at {file_path}")
        return None
    except json.JSONDecodeError:
        print(f"Error: Could not decode JSON from {file_path}")
        return None

    all_data = []

    # Iterate over each model type (e.g., 'type1', 'type2', etc.)
    for type_name, type_data in data.items():
        if FEATURE_IMPORTANCE_KEY in type_data:
            importance_data = type_data[FEATURE_IMPORTANCE_KEY]

            # Convert the dictionary of features and scores to a list of dicts
            for feature, importance in importance_data.items():
                all_data.append({
                    'Model Type': type_name,
                    'Feature': feature,
                    'Importance Score': importance
                })

    # Create the DataFrame
    df = pd.DataFrame(all_data)
    
    # Simple normalization check (optional, but good practice)
    # The sum of importance scores per type should ideally be 1.0 (or close)
    df_grouped = df.groupby('Model Type')['Importance Score'].sum()
    for model, total in df_grouped.items():
        if abs(total - 1.0) > 0.01:
             print(f"Warning: Importance scores for {model} sum to {total:.2f}. Consider normalizing if needed.")
             
    return df

def visualize_feature_importance(df):
    """
    Generates a grouped bar chart to compare feature importance across types.
    
    The top features are now selected based on the sum of their absolute
    importance scores across all model types, as requested.
    """
    
    # Calculate the total importance (sum of absolute values) across all types 
    # to determine the plot order and top features.
    feature_total_importance = df.groupby('Feature')['Importance Score'].apply(
        lambda x: x.abs().sum()
    ).sort_values(ascending=False)
    
    # Filter for the top N features based on overall total importance
    top_features = feature_total_importance.head(TOP_N_FEATURES).index.tolist()
    df_filtered = df[df['Feature'].isin(top_features)]

    # Set up the plot style
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(12, 7))

    # Create the grouped bar plot
    # The order of x-axis features is explicitly set by the total importance ranking
    ax = sns.barplot(
        x='Feature', 
        y='Importance Score', 
        hue='Model Type', 
        data=df_filtered, 
        order=top_features, # Use the calculated order
        palette='viridis',
        edgecolor=".2"
    )

    # Add labels and title
    plt.title(f'Top {TOP_N_FEATURES} Features (Ranked by Total Magnitude) Across Model Types', fontsize=16, pad=20)
    plt.xlabel('Feature Name', fontsize=12)
    plt.ylabel('Normalized Importance Score', fontsize=12)
    
    # Rotate x-axis labels for readability
    plt.xticks(rotation=45, ha='right')
    
    # Move the legend outside the plot
    plt.legend(title='Model Type', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Add context to the plot
    plt.tight_layout()
    plt.show()

# --- Main execution ---
if __name__ == '__main__':
    # 1. Load and process the data
    importance_df = load_and_process_data(JSON_FILE)

    if importance_df is not None and not importance_df.empty:
        # 2. Visualize the results
        visualize_feature_importance(importance_df)
    else:
        print("Could not generate visualization due to data loading/processing error.")