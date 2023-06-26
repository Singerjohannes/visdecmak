#!/bin/bash

# Function to unzip all zip files in a directory
unzip_files() {
    local directory=$1
    cd $directory
    for file in $directory/*.zip; do
        if [ -f "$file" ]; then
            dest=$(basename "$file" .zip)
            echo $dest
            unzip -j $file -d $dest
            rm "$file"  # Optional: Remove the zip file after extraction
        fi
    done
}

# Prompt the user for the path to the zip file
read -p "Enter the path to the zip file downloaded from the OSF fMRI data component: " zip_file_path

# Prompt the user for the destination path for extraction
read -p "Enter the path of the cloned github directory "visdecmak" for extraction: " dest_path

# Check if the zip file exists
if [ -f "$zip_file_path" ]; then
    # Unzip the file to the destination path
    unzip "$zip_file_path" -d "$dest_path"
    
    echo "Unzipping betas..."
    # Unzip all zip files within "betas" directory
    unzip_files "$dest_path/betas"
    
    echo "Unzipping deformation fields..."
    # Unzip all zip files within "deformation_fields" directory
    unzip_files "$dest_path/deformation_field"
    
    echo "Unzipping ROI files..."
    # Unzip all zip files within "roi" directory
    unzip_files "$dest_path/roi"

    echo "Unzipped successfully!"

else
    echo "Error: The specified zip file does not exist."
fi


# Prompt the user for the path to the zip file
read -p "Enter the path to the zip file downloaded from the OSF behavior component: " zip_file_path

# Check if the zip file exists
if [ -f "$zip_file_path" ]; then
    # Unzip the file to the destination path
    unzip "$zip_file_path" -d "$dest_path/behav"
    
    echo "Unzipped successfully!"

else
    echo "Error: The specified zip file does not exist."
fi