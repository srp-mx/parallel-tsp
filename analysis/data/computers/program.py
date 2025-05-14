#!/usr/bin/env python3
import os
import platform
import socket
import csv
import subprocess
import re

import psutil

# Try to import cpuinfo (py-cpuinfo)
try:
    import cpuinfo
except ImportError:
    cpuinfo = None

# Try to import NVML (for GPU info via pynvml)
try:
    import pynvml
    pynvml.nvmlInit()
except Exception:
    pynvml = None

# Optionally use PyCUDA for more detailed GPU info
try:
    import pycuda.driver as cuda
    cuda.init()
    has_pycuda = True
except Exception:
    has_pycuda = False


def get_hostname():
    return socket.gethostname()


def get_cpu_info():
    info = {}
    if cpuinfo:
        cpu_data = cpuinfo.get_cpu_info()
        info['CPU Model'] = cpu_data.get('brand_raw', 'N/A')
        info['CPU Vendor'] = cpu_data.get('vendor_id_raw', 'N/A')
        info['CPU Flags'] = " ".join(cpu_data.get('flags', []))
    else:
        model = "N/A"
        vendor = "N/A"
        flags = ""
        try:
            with open('/proc/cpuinfo') as f:
                for line in f:
                    if "model name" in line and model == "N/A":
                        model = line.split(":", 1)[1].strip()
                    if "vendor_id" in line and vendor == "N/A":
                        vendor = line.split(":", 1)[1].strip()
                    if "flags" in line and not flags:
                        flags = line.split(":", 1)[1].strip()
                        break
        except Exception:
            pass
        info['CPU Model'] = model
        info['CPU Vendor'] = vendor
        info['CPU Flags'] = flags
    return info


def get_core_info():
    info = {}
    info['Logical Cores'] = psutil.cpu_count(logical=True)
    info['Physical Cores'] = psutil.cpu_count(logical=False)
    freq = psutil.cpu_freq()
    if freq:
        info['Base Frequency (MHz)'] = freq.min
        info['Boost Frequency (MHz)'] = freq.max
    else:
        info['Base Frequency (MHz)'] = "N/A"
        info['Boost Frequency (MHz)'] = "N/A"
    return info


def get_chipset_info():
    chipset = "N/A"
    try:
        output = subprocess.check_output(
            "sudo dmidecode -t baseboard", shell=True, stderr=subprocess.DEVNULL, universal_newlines=True
        )
        m = re.search(r"Product Name:\s*(.+)", output)
        if m:
            chipset = m.group(1).strip()
    except Exception:
        chipset = "N/A"
    return chipset


def get_os_info():
    uname = platform.uname()
    distro_name = "N/A"
    try:
        with open("/etc/os-release") as f:
            for line in f:
                if line.startswith("PRETTY_NAME"):
                    distro_name = line.split("=", 1)[1].strip().strip('"')
                    break
    except Exception:
        pass
    return f"{distro_name}, Kernel: {uname.release}"


def get_glibc_version():
    libc_info = platform.libc_ver()
    if libc_info[1]:
        return libc_info[1]
    return "N/A"


def get_cache_info():
    cache_info = {}
    try:
        output = subprocess.check_output("lscpu", shell=True, universal_newlines=True)
        for line in output.splitlines():
            if "L1d cache" in line:
                cache_info['L1d Cache'] = line.split(":", 1)[1].strip()
            elif "L1i cache" in line:
                cache_info['L1i Cache'] = line.split(":", 1)[1].strip()
            elif "L2 cache" in line:
                cache_info['L2 Cache'] = line.split(":", 1)[1].strip()
            elif "L3 cache" in line:
                cache_info['L3 Cache'] = line.split(":", 1)[1].strip()
    except Exception:
        cache_info['L1d Cache'] = "N/A"
        cache_info['L1i Cache'] = "N/A"
        cache_info['L2 Cache'] = "N/A"
        cache_info['L3 Cache'] = "N/A"
    return cache_info


def get_openmp_version():
    try:
        output = subprocess.check_output("echo | gcc -fopenmp -dM -E -", shell=True, universal_newlines=True)
        m = re.search(r"#define _OPENMP\s+(\d+)", output)
        if m:
            omp_value = int(m.group(1))
            mapping = {
                200203: "2.0",
                200505: "2.5",
                200806: "3.0",
                201107: "3.1",
                201307: "4.0",
                201511: "4.5",
                202011: "5.0",
            }
            return mapping.get(omp_value, str(omp_value))
        else:
            return "N/A"
    except Exception:
        return "N/A"


def get_ram_info():
    ram_info = {}
    total_mem_bytes = psutil.virtual_memory().total
    total_mem_gb = total_mem_bytes / (1024 ** 3)
    ram_info['Total RAM (GB)'] = f"{total_mem_gb:.2f}"
    speed = "N/A"
    cas_latency = "N/A"
    try:
        output = subprocess.check_output(
            "sudo dmidecode -t memory", shell=True, stderr=subprocess.DEVNULL, universal_newlines=True
        )
        speeds = re.findall(r"Speed:\s*(\d+)\s*MHz", output)
        if speeds:
            speed = speeds[0] + " MHz"
        m = re.search(r"CAS Latency:\s*(\d+)", output)
        if m:
            cas_latency = m.group(1)
    except Exception:
        speed = "N/A"
        cas_latency = "N/A"
    ram_info['RAM Speed'] = speed
    ram_info['RAM CAS Latency'] = cas_latency
    return ram_info


def get_gpu_info():
    gpu_info = {}
    if has_pycuda:
        try:
            dev = cuda.Device(0)
            attributes = dev.get_attributes()
            gpu_info['CUDA Capable Device'] = "Yes"
            gpu_info['NVIDIA GPU Model'] = dev.name()
            total_mem = dev.total_memory()  # in bytes
            gpu_info['Global Memory'] = f"{total_mem / (1024 ** 3):.2f} GB"
            shared_mem = attributes.get(cuda.device_attribute.MAX_SHARED_MEMORY_PER_BLOCK, None)
            gpu_info['Shared Memory per Block'] = f"{shared_mem / 1024:.2f} KB" if shared_mem is not None else "N/A"
            warp_size = attributes.get(cuda.device_attribute.WARP_SIZE, None)
            gpu_info['Warp Size'] = warp_size if warp_size is not None else "N/A"
            comp_cap = dev.compute_capability()
            gpu_info['Compute Capability'] = f"{comp_cap[0]}.{comp_cap[1]}"
            mp_count = attributes.get(cuda.device_attribute.MULTIPROCESSOR_COUNT, 0)
            gpu_info['Streaming Multiprocessor Count'] = mp_count
            cc_major = comp_cap[0]
            cores_per_mp = {
                2: 32,
                3: 192,
                5: 128,
                6: 64,
                7: 64,
                8: 64,
            }.get(cc_major, "N/A")
            if isinstance(cores_per_mp, int):
                gpu_info['CUDA Core Count'] = mp_count * cores_per_mp
            else:
                gpu_info['CUDA Core Count'] = "N/A"
            gpu_info['PCIe Bus Bandwidth'] = "N/A"
        except Exception as e:
            gpu_info['CUDA Capable Device'] = f"Error retrieving GPU info: {e}"
    elif pynvml:
        try:
            device_count = pynvml.nvmlDeviceGetCount()
            if device_count > 0:
                handle = pynvml.nvmlDeviceGetHandleByIndex(0)
                gpu_info['CUDA Capable Device'] = "Yes"
                name = pynvml.nvmlDeviceGetName(handle)
                if isinstance(name, bytes):
                    name = name.decode('utf-8')
                gpu_info['NVIDIA GPU Model'] = name
                mem_info = pynvml.nvmlDeviceGetMemoryInfo(handle)
                gpu_info['Global Memory'] = f"{mem_info.total / (1024 ** 3):.2f} GB"
                gpu_info['Shared Memory per Block'] = "N/A"
                gpu_info['Warp Size'] = "N/A"
                gpu_info['Streaming Multiprocessor Count'] = "N/A"
                gpu_info['Compute Capability'] = "N/A"
                gpu_info['CUDA Core Count'] = "N/A"
                gpu_info['PCIe Bus Bandwidth'] = "N/A"
            else:
                gpu_info['CUDA Capable Device'] = "No"
        except Exception as e:
            gpu_info['CUDA Capable Device'] = f"Error retrieving GPU info: {e}"
    else:
        gpu_info['CUDA Capable Device'] = "N/A"
    try:
        cuda_ver_output = subprocess.check_output(
            "nvcc --version", shell=True, stderr=subprocess.DEVNULL, universal_newlines=True
        )
        m = re.search(r"release (\S+),", cuda_ver_output)
        if m:
            gpu_info['CUDA Version'] = m.group(1)
        else:
            gpu_info['CUDA Version'] = "N/A"
    except Exception:
        gpu_info['CUDA Version'] = "N/A"
    return gpu_info


def get_additional_cpu_extensions(flags_str):
    extensions = {}
    lower_flags = flags_str.lower() if flags_str else ""
    extensions['AVX Supported'] = "Yes" if "avx" in lower_flags else "No"
    extensions['SSE Supported'] = "Yes" if "sse" in lower_flags else "No"
    return extensions


def main():
    data = {}

    # Gather basic system info.
    data['Hostname'] = get_hostname()
    data['OS Info'] = get_os_info()
    data['glibc Version'] = get_glibc_version()

    # CPU data.
    cpu_data = get_cpu_info()
    data.update(cpu_data)
    core_data = get_core_info()
    data.update(core_data)
    data['Chipset'] = get_chipset_info()
    cache_info = get_cache_info()
    data.update(cache_info)
    data['OpenMP Version'] = get_openmp_version()

    # RAM data.
    ram_info = get_ram_info()
    data.update(ram_info)

    # GPU data.
    gpu_info = get_gpu_info()
    data.update(gpu_info)

    # Additional CPU extensions.
    if 'CPU Flags' in data:
        extra_extensions = get_additional_cpu_extensions(data['CPU Flags'])
        data.update(extra_extensions)

    # Name the CSV file after the hostname.
    hostname = data.get("Hostname", "benchmark_info")
    csv_file = f"{hostname}.csv"

    try:
        with open(csv_file, mode="w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["Parameter", "Value"])
            for key, value in data.items():
                writer.writerow([key, value])
        print(f"Benchmark information written to {csv_file}")
    except Exception as e:
        print(f"Failed to write CSV file: {e}")


if __name__ == "__main__":
    main()

