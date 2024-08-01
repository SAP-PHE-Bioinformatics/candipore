import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import concurrent.futures

def process_bedgraph(file_path):
    # Read bedGraph file into a DataFrame
    df = pd.read_csv(file_path, sep='\t', header=None, names=['chrom', 'start', 'end', 'coverage'])
    
    # Calculate the length of each interval
    df['length'] = df['end'] - df['start']
    
    # Calculate total genome length covered by bedGraph
    total_length = df['length'].sum()
    
    # Calculate coverage statistics
    zero_coverage = df[df['coverage'] == 0]['length'].sum() / total_length * 100
    less_than_10_reads = df[(df['coverage'] > 0) & (df['coverage'] < 10)]['length'].sum() / total_length * 100
    more_than_10_reads = df[df['coverage'] >= 10]['length'].sum() / total_length * 100
    
    return {
        'file_path': file_path,
        'Sample': file_path.split('/')[-1],
        'Zero Coverage (%)': zero_coverage,
        'Coverage < 10 Reads (%)': less_than_10_reads,
        'Coverage >= 10 Reads (%)': more_than_10_reads,
        'coordinates': df['start'].tolist(),
        'coverages': df['coverage'].tolist()
    }

def plot_coverage(stats):
    coordinates = stats['coordinates']
    coverages = stats['coverages']
    
    # Plot the coverage data
    plt.figure(figsize=(10, 6))
    plt.plot(coordinates, coverages, linestyle='none', marker=',')
    plt.xlabel('Coordinate')
    plt.ylabel('Number of Reads')
    plt.title(f'Coverage Plot for {stats["Sample"]}')
    plt.show()
    plt.savefig(f'coverage_plot.png')

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Process and plot coverage statistics from a bedGraph file.')
    parser.add_argument('--bedgraph', nargs='+', type=str, help='Path to one or more bedGraph files')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Process each bedGraph file in parallel
    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = list(executor.map(process_bedgraph, args.bedgraph))
    
    # Create a DataFrame from the statistics
    stats_list = [result for result in results]
    stats_df = pd.DataFrame(stats_list).drop(columns=['coordinates', 'coverages'])
    
    # Display the DataFrame
    print(stats_df)
    
    # Optionally save the DataFrame to a CSV file
    stats_df.to_csv('coverage_stats.csv', index=False)
    
    # Plot coverage for each bedGraph file (sequentially in the main thread)
    for stats in results:
        plot_coverage(stats)

if __name__ == '__main__':
    main()
