


# Function name: obtainRawData
# Author: Alberto Poncelas

# This function downloads .CEL files from 
# Gene Expression Omnibus (see http://www.ncbi.nlm.nih.gov/sites/GDSbrowser) and 
# extract the .CEL files into a folder. (for just downloading the 
# compressed .tar file, use "getGEOSuppFiles" command)

# The function parameters are:
# ----accession_number:   The accession number of the data set 
#				(the code of "series" column in http://www.ncbi.nlm.nih.gov/sites/GDSbrowser)
# ----folder_name:		The name of the folder where .CEL files will be stored

# Note: This function has been tested only with some raw data. Is not garanteed that to work with every data.
# This function creates and deletes temporal folders and files, so it is recommended to use in project's folder



obtainRawData<-function(accession_number,folder="data"){

	library(GEOquery)
	library(R.utils)
	
	
	#Define the path of a temporal folder to extract files
	temp_folder=paste(getwd(),"/temp",sep="")
	
	#Define and create the path to save .CEL files
	data_folder=paste(getwd(),folder,sep="/")
	dir.create(data_folder)


	#Download .tar with files from GEO
	getGEOSuppFiles(GEO=accession_number)


	#Obtain the path of rae data and unzip it in temporal folder
	raw_file_name=paste(accession_number,"_RAW.tar",sep="")
	raw_file_folder=paste(getwd(),"/",accession_number,"/",raw_file_name,sep="")
	untar(raw_file_folder, exdir=temp_folder)

	#Get the names of .cel.gz compressed files
	cel_zip_files <- list.files(temp_folder, pattern = "cel.gz$",ignore.case=TRUE)

	#Unzip the files that contains .CEL files
	sapply(paste(temp_folder, cel_zip_files , sep="/"), gunzip, remove=TRUE)


	#Get the names of .CEL files
	cel_files <- list.files(temp_folder, pattern = ".cel$",ignore.case=TRUE)


	#Copy those files from temporal folder to data folder
	files_from<-paste(temp_folder, cel_files, sep="/")
	files_to<-paste(data_folder, cel_files, sep="/")
	file.copy(from=files_from,to=files_to)


	#Delete temporal folder
	unlink(temp_folder, recursive = TRUE)

	#Delete downloaded .tar with raw data
	raw_file_folder=paste(getwd(),"/",accession_number,sep="")
	unlink(raw_file_folder, recursive = TRUE)

}

