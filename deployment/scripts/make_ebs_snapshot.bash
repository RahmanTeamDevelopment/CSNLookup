#!/bin/bash

# Create an EBS volume large enough to store the 2 reference files, indexes, and possibly other data
aws ec2 create-volume --size 50 --region eu-west-2 --availability-zone eu-west-2a > ebs_volume.json

# Create an EC2 instance to attach the volume to (the ami is one I pulled from the AWS console. It runs a small instance with Amazon Linux)
aws ec2 run-instances --image-id ami-2d09eb5a --count 1 --instance-type t2.small

# Attach volume to EC2 instance
aws ec2 attach-volume --volume-id $EBS_VOLUME_ID --instance-id $EC2_INSTANCE_ID --device /dev/sdf

# Mount the EBS volume as a file-system. These commands must be run on the EC2 instance (remotely with AWS CLI or Fabric?)
sudo mkfs -t ext3 /dev/sdh
sudo mkdir /media/volume1
sudo mount /dev/sdh /media/volume1

# Copy reference files to EC2 instance. These should also be run on the EC2 instance to take advantage of the fast connections between S3 and EBS volumes
# in the Amazon data centre
aws s3 cp s3://andy-temp-storage/hg38.fa /media/volume1/
aws s3 cp s3://andy-temp-storage/hg38.fa.fai /media/volume1/
aws s3 cp s3://andy-temp-storage/human_g1k_v37.fa /media/volume1/
aws s3 cp s3://andy-temp-storage/human_g1k_v37.fa.fai /media/volume1/

# Unmount EBS volume from EC2 instance in preparation for creating snapshot (run on EC2 instance)
sudo umount /media/volume1

# Create snapshot
aws ec2 create-snapshot --volume-id $EBS_VOLUME_ID --description "Snap-shot of EBS volume storing 37 and 37 reference files"

# Delete the EBS volume. We don't need it once we have taken the snapshot
aws ec2 delete-volume $EBS_VOLUME_ID

# Terminate the EC2 instance
aws ec2 terminate-instances --instance-ids $EC2_INSTANCE_ID
