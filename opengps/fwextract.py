#!/usr/bin/python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Cut the MT6755 GPS-core firmware out of the vendor's libmnl.so (like b43-fwcutter: nothing
proprietary is redistributed; the image is taken from the device's own vendor partition).

usage: fwextract.py [path/to/libmnl.so] [out=/lib/firmware/mediatek/mt6755_gps.bin]
       (default: the first existing of /vendor/lib/libmnl.so, /system/vendor/lib/libmnl.so,
        /mnt/vendor/lib/libmnl.so - mount the device's vendor partition there)

libmnl.so carries several images that start with the same bytes (one per chip); the one stock
mnld downloads on this tablet (44 x FE0C chunks, 61440 bytes) is identified by its SHA-256.
"""
import hashlib, os, sys

KNOWN = {  # sha256 -> (size, description)
    "ac35cea91bfa03bd8158411576d01d51644c699603855a8cc1db0bf2e4b2a332": (61440, "MT6755 (libmnl.so, mnld 2008-12-31 build)"),
}
HEAD = bytes.fromhex("fb417bf3fb41c982fb414681fb412681")


def find(lib):
    i = -1
    while True:
        i = lib.find(HEAD, i + 1)
        if i < 0:
            return None
        for sha, (size, desc) in KNOWN.items():
            img = lib[i:i + size]
            if hashlib.sha256(img).hexdigest() == sha:
                return i, img, desc


def main():
    cands = ["/vendor/lib/libmnl.so", "/system/vendor/lib/libmnl.so", "/mnt/vendor/lib/libmnl.so"]
    src = sys.argv[1] if len(sys.argv) > 1 else next((c for c in cands if os.path.exists(c)), None)
    if not src:
        sys.exit("libmnl.so not found; pass the path of the device's vendor/lib/libmnl.so")
    dst = sys.argv[2] if len(sys.argv) > 2 else "/lib/firmware/mediatek/mt6755_gps.bin"
    r = find(open(src, "rb").read())
    if not r:
        sys.exit("no known GPS firmware image in %s" % src)
    off, img, desc = r
    os.makedirs(os.path.dirname(dst) or ".", exist_ok=True)
    open(dst, "wb").write(img)
    print("%s: %d bytes at offset %d (%s) -> %s" % (src, len(img), off, desc, dst))


if __name__ == "__main__":
    main()
