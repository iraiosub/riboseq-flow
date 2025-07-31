#!/usr/bin/env nextflow

def getAdapterChannel(adapterFile, noAdapterValue) {
    return adapterFile ? Channel.fromPath(adapterFile, checkIfExists: true).first()
                       : Channel.value(noAdapterValue).collectFile{ it -> ["${it}", it] }.first()
}
