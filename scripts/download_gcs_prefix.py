from __future__ import annotations

import argparse
import json
import os
import urllib.parse
import urllib.request
from pathlib import Path


TOKEN_URL = "https://oauth2.googleapis.com/token"
STORAGE_LIST_URL = "https://storage.googleapis.com/storage/v1/b/{bucket}/o"
STORAGE_DOWNLOAD_URL = "https://storage.googleapis.com/download/storage/v1/b/{bucket}/o/{object_name}?alt=media"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Download all objects under a GCS prefix using ADC credentials.")
    parser.add_argument("--bucket", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument(
        "--credentials",
        default=os.path.expanduser("~/.config/gcloud/application_default_credentials.json"),
        help="Path to application_default_credentials.json",
    )
    return parser


def fetch_access_token(credentials_path: str) -> str:
    with open(credentials_path, "r", encoding="utf-8") as handle:
        creds = json.load(handle)
    payload = urllib.parse.urlencode(
        {
            "client_id": creds["client_id"],
            "client_secret": creds["client_secret"],
            "refresh_token": creds["refresh_token"],
            "grant_type": "refresh_token",
        }
    ).encode("utf-8")
    request = urllib.request.Request(TOKEN_URL, data=payload, method="POST")
    request.add_header("Content-Type", "application/x-www-form-urlencoded")
    with urllib.request.urlopen(request, timeout=60) as response:
        token_payload = json.loads(response.read().decode("utf-8"))
    return token_payload["access_token"]


def list_objects(bucket: str, prefix: str, access_token: str) -> list[dict[str, str]]:
    items: list[dict[str, str]] = []
    page_token: str | None = None
    while True:
        params = {"prefix": prefix}
        if page_token is not None:
            params["pageToken"] = page_token
        url = STORAGE_LIST_URL.format(bucket=bucket) + "?" + urllib.parse.urlencode(params)
        request = urllib.request.Request(url)
        request.add_header("Authorization", f"Bearer {access_token}")
        with urllib.request.urlopen(request, timeout=120) as response:
            payload = json.loads(response.read().decode("utf-8"))
        items.extend(payload.get("items", []))
        page_token = payload.get("nextPageToken")
        if not page_token:
            break
    return items


def download_object(bucket: str, object_name: str, output_dir: Path, access_token: str) -> Path:
    encoded_name = urllib.parse.quote(object_name, safe="")
    url = STORAGE_DOWNLOAD_URL.format(bucket=bucket, object_name=encoded_name)
    request = urllib.request.Request(url)
    request.add_header("Authorization", f"Bearer {access_token}")
    destination = output_dir / object_name
    destination.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(request, timeout=600) as response, destination.open("wb") as handle:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            handle.write(chunk)
    return destination


def main() -> None:
    args = build_parser().parse_args()
    output_dir = Path(args.output_dir).resolve()
    access_token = fetch_access_token(args.credentials)
    items = list_objects(args.bucket, args.prefix, access_token)
    print(f"Found {len(items)} objects under gs://{args.bucket}/{args.prefix}")
    for item in items:
        object_name = item["name"]
        size = int(item.get("size", 0))
        destination = download_object(args.bucket, object_name, output_dir, access_token)
        print(f"Downloaded {object_name} -> {destination} ({size} bytes)")


if __name__ == "__main__":
    main()
