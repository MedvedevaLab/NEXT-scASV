#!/usr/bin/env python3

"""
Utility functions for parsing meta.json files in the ASE pipeline
"""

import json
import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any


class MetaParser:
    @staticmethod
    def parse_meta_json(meta_json_path: str) -> Dict:
        """
        Parse meta.json file and return the parsed JSON object. Repalce "-" in sample_id with "_"
        
        Args:
            meta_json_path: Path to the meta.json file
            
        Returns:
            Parsed JSON object
        """
        with open(meta_json_path, 'r') as f:
            metadata = json.load(f)
            for sample_id in metadata.keys():
                metadata[sample_id]['groups'] = [group.replace('-', '_') for group in metadata[sample_id]['groups']]
            return metadata

    @staticmethod
    def get_sample_ids(metadata: Dict) -> List[str]:
        """
        Get list of sample IDs from meta.json
        
        Args:
            metadata: Parsed meta.json object
            
        Returns:
            List of sample IDs
        """
        return list(metadata.keys())

    @staticmethod
    def create_entries_list(metadata: Dict) -> List[Dict[str, str]]:
        """
        Create entries list from meta.json for processing
        
        Args:
            metadata: Parsed meta.json object
            
        Returns:
            List of entries with sample_id and group information
        """
        entries = []
        for sample_id, sample_info in metadata.items():
            for group in sample_info['groups']:
                entries.append({
                    'sample_id': sample_id,
                    'group': group
                })
        return entries

    @staticmethod
    def create_sample_name(sample_id: str, group: str) -> str:
        """
        Create sample name from sample_id and group
        
        Args:
            sample_id: Sample ID
            group: Group name
            
        Returns:
            Formatted sample name
        """
        return f"{sample_id}-{group}"

    @staticmethod
    def create_bam_entries(
        metadata: Dict,
        bam_dir: str,
        suffix: str = "_dedup_with_tag.bam"
    ) -> List[Tuple[str, str, str, str]]:
        """
        Create BAM entries from metadata
        
        Args:
            metadata: Parsed meta.json object
            bam_dir: Directory containing BAM files
            suffix: Suffix to append to sample name for BAM file
            
        Returns:
            List of BAM entries with sample information
        """
        bam_entries = []
        for sample_id, sample_info in metadata.items():
            for group in sample_info['groups']:
                sample_name = MetaParser.create_sample_name(sample_id, group)
                bam_path = Path(bam_dir) / f"{sample_name}{suffix}"
                bai_path = Path(str(bam_path) + '.bai')
                
                if not bam_path.exists():
                    print(f"WARNING: BAM file not found: {bam_path}", file=sys.stderr)
                    continue
                
                if not bai_path.exists():
                    print(f"WARNING: Index file not found for {bam_path}", file=sys.stderr)
                    continue
                
                bam_entries.append([sample_name, sample_id, str(bam_path), str(bai_path)])
        return bam_entries

    @staticmethod
    def create_bam_bed_entries(
        metadata: Dict,
        bam_dir: str,
        bed_dir: str,
        suffix: str = "_marked_filtered.bam"
    ) -> List[Tuple[str, str, str, str, str, str]]:
        """
        Create BAM entries with BED files from metadata
        
        Args:
            metadata: Parsed meta.json object
            bam_dir: Directory containing BAM files
            bed_dir: Directory containing BED files
            suffix: Suffix to append to sample name for BAM file
            
        Returns:
            List of entries with sample information and associated files
        """
        entries = []
        for sample_id, sample_info in metadata.items():
            for group in sample_info['groups']:
                sample_name = MetaParser.create_sample_name(sample_id, group)
                bam_path = Path(bam_dir) / f"{sample_name}{suffix}"
                bai_path = Path(str(bam_path) + '.bai')
                bed_path = Path(bed_dir) / f"{sample_id}.bed.gz"
                tbi_path = Path(str(bed_path) + '.tbi')
                
                if not bam_path.exists():
                    print(f"WARNING: BAM file not found: {bam_path}", file=sys.stderr)
                    continue
                
                if not bai_path.exists():
                    print(f"WARNING: Index file not found for {bam_path}", file=sys.stderr)
                    continue

                if not bed_path.exists() or not tbi_path.exists():
                    print(f"WARNING: BED or TBI file not found for {sample_id}: {bed_path}", file=sys.stderr)
                    continue
                
                entries.append([sample_name, sample_id, str(bam_path), str(bai_path), str(bed_path), str(tbi_path)])
        return entries

    @staticmethod
    def create_sample_map_entries(metadata: Dict) -> List[Tuple[str, str, str]]:
        """
        Create sample map entries for recode_vcf process from metadata
        
        Args:
            metadata: Parsed meta.json object
            
        Returns:
            List of sample map entries with sample_name, sample_id, and expected counts file name
        """
        entries = []
        for sample_id, sample_info in metadata.items():
            for group in sample_info['groups']:
                sample_name = MetaParser.create_sample_name(sample_id, group)
                counts_file_name = f"{sample_name}.counts.bed.gz"
                entries.append([sample_name, sample_id, counts_file_name])
        return entries

    @staticmethod
    def get_unique_groups(metadata: Dict) -> List[str]:
        """
        Extract unique groups from metadata
        
        Args:
            metadata: Parsed meta.json object
            
        Returns:
            List of unique group names
        """
        groups = set()
        for sample_id, sample_info in metadata.items():
            groups.update(sample_info['groups'])
        return sorted(list(groups))

    @staticmethod
    def create_ase_file_names(metadata: Dict, counts_dir: str) -> List[str]:
        """
        Create list of count file paths for ASE analysis from metadata
        
        Args:
            metadata: Parsed meta.json object
            counts_dir: Directory containing count files
            
        Returns:
            List of count file paths
        """
        file_names = []
        for sample_id, sample_info in metadata.items():
            for group in sample_info['groups']:
                sample_name = MetaParser.create_sample_name(sample_id, group)
                count_file_path = Path(counts_dir) / f"{sample_name}.counts.bed.gz"
                if count_file_path.exists():
                    file_names.append(str(count_file_path))
                else:
                    print(f"WARNING: Count file not found: {count_file_path}", file=sys.stderr)
        return file_names

    @staticmethod
    def filter_vcf_files_by_group(metadata: Dict, vcf_files: List[str], group: str) -> List[str]:
        """
        Filter VCF files by group based on metadata
        
        Args:
            metadata: Parsed meta.json object
            vcf_files: List of VCF file paths
            group: Group name to filter by
            
        Returns:
            List of VCF file paths that belong to the specified group
        """
        filtered_files = []
        for vcf_file in vcf_files:
            # Extract sample name from VCF filename (remove .vcf.gz extension)
            vcf_basename = Path(vcf_file).name.replace('.vcf.gz', '').replace('.vcf', '')
            # Check if this VCF file matches the group pattern (sample_id-group)
            if vcf_basename.endswith(f"-{group}"):
                filtered_files.append(vcf_file)
        return filtered_files


def main():
    if len(sys.argv) < 2:
        print("Usage: meta_parser.py <command> [args...]", file=sys.stderr)
        sys.exit(1)

    command = sys.argv[1]
    
    if command == "parse_meta_json":
        if len(sys.argv) != 3:
            print("Usage: meta_parser.py parse_meta_json <meta_json_path>", file=sys.stderr)
            sys.exit(1)
        metadata = MetaParser.parse_meta_json(sys.argv[2])
        print(json.dumps(metadata, indent=None, separators=(',', ':')))
        
    elif command == "get_sample_ids":
        if len(sys.argv) != 3:
            print("Usage: meta_parser.py get_sample_ids <meta_json_path>", file=sys.stderr)
            sys.exit(1)
        metadata = MetaParser.parse_meta_json(sys.argv[2])
        sample_ids = MetaParser.get_sample_ids(metadata)
        print(json.dumps(sample_ids, indent=None, separators=(',', ':')))
        
    elif command == "create_entries_list":
        if len(sys.argv) != 3:
            print("Usage: meta_parser.py create_entries_list <meta_json_path>", file=sys.stderr)
            sys.exit(1)
        metadata = MetaParser.parse_meta_json(sys.argv[2])
        entries = MetaParser.create_entries_list(metadata)
        # Output each entry on a new line for Nextflow to process
        for entry in entries:
            print(f"{entry['sample_id']}\t{entry['group']}")
        
    elif command == "create_sample_name":
        if len(sys.argv) != 4:
            print("Usage: meta_parser.py create_sample_name <sample_id> <group>", file=sys.stderr)
            sys.exit(1)
        sample_name = MetaParser.create_sample_name(sys.argv[2], sys.argv[3])
        print(sample_name)
        
    elif command == "create_bam_entries":
        if len(sys.argv) != 4:
            print("Usage: meta_parser.py create_bam_entries <meta_json_path> <bam_dir>", file=sys.stderr)
            sys.exit(1)
        metadata = MetaParser.parse_meta_json(sys.argv[2])
        entries = MetaParser.create_bam_entries(metadata, sys.argv[3])
        # Output each entry on a new line with tab separation
        for entry in entries:
            print("\t".join(entry))
        
    elif command == "create_bam_bed_entries":
        if len(sys.argv) != 5:
            print("Usage: meta_parser.py create_bam_bed_entries <meta_json_path> <bam_dir> <bed_dir>", file=sys.stderr)
            sys.exit(1)
        metadata = MetaParser.parse_meta_json(sys.argv[2])
        entries = MetaParser.create_bam_bed_entries(metadata, sys.argv[3], sys.argv[4])
        # Output as JSON for NextFlow to parse
        print(json.dumps(entries, indent=None, separators=(',', ':')))
        
    elif command == "create_sample_map_entries":
        if len(sys.argv) != 3:
            print("Usage: meta_parser.py create_sample_map_entries <meta_json_path>", file=sys.stderr)
            sys.exit(1)
        metadata = MetaParser.parse_meta_json(sys.argv[2])
        entries = MetaParser.create_sample_map_entries(metadata)
        # Output each entry on a new line with tab separation
        for entry in entries:
            print("\t".join(entry))
        
    elif command == "get_unique_groups":
        if len(sys.argv) != 3:
            print("Usage: meta_parser.py get_unique_groups <meta_json_path>", file=sys.stderr)
            sys.exit(1)
        metadata = MetaParser.parse_meta_json(sys.argv[2])
        groups = MetaParser.get_unique_groups(metadata)
        print(json.dumps(groups, indent=None, separators=(',', ':')))
        
    elif command == "create_ase_file_names":
        if len(sys.argv) != 4:
            print("Usage: meta_parser.py create_ase_file_names <meta_json_path> <counts_dir>", file=sys.stderr)
            sys.exit(1)
        metadata = MetaParser.parse_meta_json(sys.argv[2])
        file_names = MetaParser.create_ase_file_names(metadata, sys.argv[3])
        # Output each file name on a new line
        for file_name in file_names:
            print(file_name)
        
    elif command == "get_group_ids":
        if len(sys.argv) != 3:
            print("Usage: meta_parser.py get_group_ids <meta_json_path>", file=sys.stderr)
            sys.exit(1)
        metadata = MetaParser.parse_meta_json(sys.argv[2])
        groups = MetaParser.get_unique_groups(metadata)
        # Output each group ID on a new line for bash loop processing
        for group in groups:
            print(group)
        
    elif command == "filter_vcf_files_by_group":
        if len(sys.argv) < 5:
            print("Usage: meta_parser.py filter_vcf_files_by_group <meta_json_path> <group> <vcf_file1> [vcf_file2] ...", file=sys.stderr)
            sys.exit(1)
        metadata = MetaParser.parse_meta_json(sys.argv[2])
        group = sys.argv[3]
        vcf_files = sys.argv[4:]
        filtered_files = MetaParser.filter_vcf_files_by_group(metadata, vcf_files, group)
        # Output each file name on a new line
        for file_name in filtered_files:
            print(file_name)
        
    else:
        print(f"Unknown command: {command}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main() 