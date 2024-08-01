import polars as pl
import matplotlib.pyplot as plt
import argparse
import concurrent.futures

def process_bedgraph(file_path):
    # Read bedGraph file into a DataFrame using Polars
    df = pl.read_csv(file_path, separator='\t', has_header=False, new_columns=['chrom', 'start', 'end', 'coverage'])
    
    # Calculate the length of each interval
    df = df.with_columns((df['end'] - df['start']).alias('length'))
    
    # Calculate total genome length covered by bedGraph
    total_length = df['length'].sum()
    
    # Calculate coverage statistics
    zero_coverage = df.filter(pl.col('coverage') == 0)['length'].sum() / total_length * 100
    less_than_10_reads = df.filter((pl.col('coverage') > 0) & (pl.col('coverage') < 10))['length'].sum() / total_length * 100
    more_than_10_reads = df.filter(pl.col('coverage') >= 10)['length'].sum() / total_length * 100
    
    # Flatten the coverage data for plotting
    coordinates = []
    coverages = []
    for row in df.iter_rows(named=True):
        coordinates.extend(range(row['start'], row['end']))
        coverages.extend([row['coverage']] * (row['end'] - row['start']))
    
    return {
        'file_path': file_path,
        'Sample': file_path.split('/')[-1],
        'Zero Coverage (%)': zero_coverage,
        'Coverage < 10 Reads (%)': less_than_10_reads,
        'Coverage >= 10 Reads (%)': more_than_10_reads,
        'coordinates': coordinates,
        'coverages': coverages
    }

def plot_coverage(stats, combined=False):
    coordinates = stats['coordinates']
    coverages = stats['coverages']
    
    plt.plot(coordinates, coverages, linestyle='none', marker=',', label=stats['Sample'])
    if not combined:
        plt.xlabel('Coordinate')
        plt.ylabel('Number of Reads')
        plt.title(f'Coverage Plot for {stats["Sample"]}')
        plt.savefig(f'coverage_plot_{stats["Sample"]}.png')

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Process and plot coverage statistics from bedGraph files.')
    parser.add_argument('--bedgraph', nargs='+', type=str, help='Path to one or more bedGraph files')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Process each bedGraph file in parallel
    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = list(executor.map(process_bedgraph, args.bedgraph))
    
    # Create a DataFrame from the statistics
    stats_list = [result for result in results]
    stats_df = pl.DataFrame(stats_list).select(['Sample', 'Zero Coverage (%)', 'Coverage < 10 Reads (%)', 'Coverage >= 10 Reads (%)'])
    
    # Display the DataFrame
    print(stats_df)
    
    # Optionally save the DataFrame to a CSV file
    stats_df.write_csv('coverage_stats.csv')
    
    # Plot coverage for each bedGraph file individually
    for stats in results:
        plot_coverage(stats)
    
    # Plot combined coverage
    plt.figure(figsize=(10, 6))
    for stats in results:
        plot_coverage(stats, combined=True)
    plt.xlabel('Coordinate')
    plt.ylabel('Number of Reads')
    plt.title('Combined Coverage Plot')
    plt.legend()
    plt.savefig('combined_coverage_plot.png')

if __name__ == '__main__':
    main()
