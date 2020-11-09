from sklearn import metrics
import numpy as np
import math

class metricor:
    def __init__(self, a = 1, probability = True, bias = 'flat', ):
        self.a = a
        self.probability = probability
        self.bias = bias 
    
    def detect_model(self, model, label, contamination = 0.1, window = 100, is_A = False):
        score = self.scale(model.decision_scores_, contamination = contamination)
        if is_A is False:
            scoreX = np.zeros(len(score)+window)
            scoreX[math.ceil(window/2): len(score)+window - math.floor(window/2)] = score 
        else:
            scoreX = score
        L = self.metric(label, scoreX)
        return L
        
    def metric(self, label, score):
        '''input prediction and labels, output list of scores on different metrics
            contamination of the data
        '''
        
        #area under curve
        if np.sum(label) == 0 or np.isnan(score).any() or score is None:
            L = [None, None, None, None, None, None]
            print('Label must have groud truth value for calculating AUC score and score \
                must not be none. ')
        else:
            auc = metrics.roc_auc_score(label, score)
            
            #precision, recall, F2, supp
            preds = np.zeros(score.shape[0])
            preds[score >= 0.5] = 1
            Precision, Recall, F2, Support = metrics.precision_recall_fscore_support(label, preds)
            precision = Precision[1]
            recall = Recall[1]
            f2 = F2[1]

            #range anomaly 
            a = self.a
            Rrecall = self.range_recall(label, score, a)

            Rprecision = self.range_recall(label, score, a=0)

            L = [auc, precision, recall, f2, Rrecall, Rprecision]

            return L
    def range_precision(self, labels, preds):
        probability = self.probability
        a = self.a
        if probability == True:
            p = self.labels_conv(preds)
            p0 = self.labels_conv_binary(preds)
            range_detect = self.range_convers(p0)       
        
    def range_recall(self, labels, preds, a):
        '''
        labels = 'lengh of dataset, binary values'
        '''
        probability = self.probability
        if probability == True:
            p = self.labels_conv(preds)
            p0 = self.labels_conv_binary(preds)
            range_detect = self.range_convers(p0)

        range_label = self.range_convers(labels)
        max_score = np.sum(labels)
        length = len(range_label)
        ## set a ranged anomoly labels 
        exist = self.existence_reward(range_label, p, max_score)


        Overlap_reward = 0
        for i in range_label:
            Overlap_reward += self.w(i, p) * self.Cardinality_factor(i, range_detect)


        score = a * exist + (1-a)*Overlap_reward
        if length != 0:
            return score/length
        else:
            return 0
        
    def labels_conv(self, preds):
        '''return indices of predicted anomaly
        '''

        p = np.zeros(len(preds))
        index = np.where(preds >= 0.5)
        return index[0]
    
    def labels_conv_binary(self, preds):
        '''return predicted label
        '''
        p = np.zeros(len(preds))
        index = np.where(preds >= 0.5)
        p[index[0]] = 1
        return p 
        
    def range_convers(self, label):
        '''
        input: arrays of binary values 
        output: list of ordered pair [[a0,b0], [a1,b1]... ] of the inputs
        '''
        L = []
        i = 0
        j = 0 
        while j < len(label):
            while label[i] != 1:
                i+=1
                if i >= len(label):
                    break
            j = i+1
            if j >= len(label):
                L.append((i,j-1))

                break
            while label[j] != 0:
                j+=1
                if j >= len(label):
                    L.append((i,j-1))
                    break
            if j >= len(label):
                break
            L.append((i, j-1))
            i = j
        return L[:-1]
    
    def existence_reward(self, labels, preds, max_score):
        '''
        labels: list of ordered pair 
        preds predicted data
        '''
        tot_score = len(labels)
        score = 0
        for i in labels:
            if np.sum(np.multiply(preds <= i[1], preds >= i[0])) > 0:
                
                score += 1
        return score

    def w(self, AnomalyRange, p):
        MyValue = 0
        MaxValue = 0
        start = AnomalyRange[0]
        AnomalyLength = AnomalyRange[1] - AnomalyRange[0] + 1
        for i in range(start, start +AnomalyLength):
            bi = self.b(i, AnomalyLength)
            MaxValue +=  bi
            if i in p:
                MyValue += bi
        return MyValue/MaxValue

    def Cardinality_factor(self, Anomolyrange, Prange):
        score = 0 
        start = Anomolyrange[0]
        end = Anomolyrange[1]
        for i in Prange:
            if i[0] >= start and i[0] <= end:
                score +=1 
            elif start >= i[0] and start <= i[1]:
                score += 1
            elif end >= i[0] and end <= i[1]:
                score += 1
            elif start >= i[0] and end <= i[1]:
                score += 1
        if score == 0:
            return 0
        else:
            return 1/score
        
    def b(self, i, length):
        bias = self.bias 
        if bias == 'flat':
            return 1
        elif bias == 'front-end bias':
            return length - i + 1
        elif bias == 'back-end bias':
            return i
        else:
            if i <= length/2:
                return i
            else:
                return length - i + 1

        
    def scale(self, score, contamination):
        threshold = np.percentile(score, 100 * (1-contamination))
        minimum = np.min(score)
        if threshold != minimum:
            score = 0.5 * (score - minimum)/(threshold - minimum)
        else:
            try:
                min2 = np.min(score[score > minimum])
                score[score > minimum] = 1
                score[score >= 1] = 1
            except:
                score = None
        #pred = score >= 0.5
        return score
        
        