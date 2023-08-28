# Use the official Python image as the base image
FROM python:3.8

# Set the working directory within the container
WORKDIR /app

# Copy the Python script and input files to the container
COPY 1-txt-to-fasta.py /app/
COPY my_sequences /app/my_sequences/

# Install the required packages
RUN pip install biopython

# Run the Python script when the container starts
CMD ["python", "1-txt-to-fasta.py"]