#!/usr/bin/sh

################################################
#File Name: Run.combine.call.orfs.20250328.sh
#Author: LiTing
#Mail: liting@stu.pku.edu.cn
#Created Time: Fri 28 Mar 2025 08:32:13 PM CST
################################################

set -eo pipefail

bash 2.2a.Run.filter_and_merge.sh

nohup bash 2.2b.Run.combine.call.orfs.20250331.sh &> log/2.2b.Run.combine.call.orfs.20250331.log &