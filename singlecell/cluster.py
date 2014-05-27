from cluster_helper import cluster as ipc

def get_cluster_view(args):
    return ipc.cluster_view(args.scheduler.lower(), args.queue,
                          args.num_jobs, args.cores_per_job,
                          extra_params={"resources": args.resources,
                                        "mem": args.memory_per_job,
                                        "tag": "singlecell",
                                        "run_local": args.local})

def wait_until_complete(jobs):
    return [j.get() for j in jobs]
