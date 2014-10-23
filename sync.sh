#!/bin/bash
rsync -trz --progress --delete --exclude='.git' ~/Dropbox/Projects/chipnexus_ec2/ chipnexus:/data/analysis_code/
