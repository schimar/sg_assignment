#! /bin/python3

# this script reads two vcf files and calculates the length of the intersection (TP) and number of loci that occur in only either of the two vcf files, respectively (FN & FP). Additionally, the sensitivity, precision and specifity are being calculated from the above values (TP, FP, & FN). The values for TP, FN, FP & TN are output as a confusion matrix, with the latter three values (sensitivity, precision & specificity) printed below.

# Usage: ./getConfusionMatrix.py <actual>.vcf <predicted>.vcf

from sys import argv
from tabulate import tabulate



def getConfMatrix(actualSet, predictedSet, bpActual):
    # intersection of vars (true positive)
    tp = len(actualSet & predictedSet)
    # vars in set1 but not in set2 (and vice versa; false negative & false positive, respectively)
    fn = len(actualSet - predictedSet)
    fp = len(predictedSet - actualSet)
    # true negative
    tn = bpActual - (tp + fp + fn)
    # derivations from confusion matrix:
    sensitivity = 100 * tp / (tp + fn)
    precision = 100 * tp / (tp + fp)
    specificity = 100 * tn / (tn + fp)
    #
    return tabulate([["actual yes", tp, fn], ["actual no", fp, tn]], ["", "predicted yes", "predicted no"], tablefmt="grid"), sensitivity, precision, specificity






if __name__ == "__main__":

    with open(argv[1], 'rb') as actual:
        actualSet = set()
        for line in actual:
            line = line.decode('utf-8')
            if line[0:1] == '#':
                continue
            else:
                linels = line.split('\t')
                bp, ref, alt = linels[1], linels[3], linels[4]
                bprefalt = ':'.join([bp, ref, alt])
                actualSet.add(bprefalt)
    actual.close()

    with open(argv[2], 'rb') as predicted:
        predictedSet = set()
        for line in predicted:
            line = line.decode('utf-8')
            if line[0:1] == '#':
                continue
            else:
                linels = line.split('\t')
                bp, ref, alt = linels[1], linels[3], linels[4]
                bprefalt = ':'.join([bp, ref, alt])
                predictedSet.add(bprefalt)
    predicted.close()




    # length (in bp) for true negatives
    # Here, I take the position of last variant - the position of first variant in ground_truth.vcf as the total length. From the task: "Any position not included in this file is assumed to be negative (i.e., no variants present)". I could either take the full length of chromosome 19, but this seems counterintuitive to me, as the alignment doesn't span across the whole chromosome 19.

    bpGroundTruth = 14999021 - 259464
    ##  14739557

    confusionMatrix, sensitivity, precision, specificity = getConfMatrix(actualSet, predictedSet, bpGroundTruth)
    print(confusionMatrix, '\n')
    print('Sensitivity = %d' % sensitivity)
    print('Precision = %d' % precision)
    print('Specificity = %d' % specificity)


