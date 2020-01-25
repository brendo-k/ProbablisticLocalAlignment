import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class Alignment {
    public String genomeAlignment;
    public String queryAlignment;
    public double probability;
    //inclusive range
    public int[] queryRangeAlignment;
    public int[] genomeRangeAignment;

    public Alignment(String sequence, int[] genomeRange, String query, int[] queryRange, Double[] prob){
        ArrayList<ArrayList<Pointer>> dpTable = new ArrayList<>();
        setTableSize(genomeRange[1] - genomeRange[0] + 2, queryRange[1] - queryRange[0] + 2, dpTable);

        Pointer biggestValue = SW(dpTable, sequence, query, genomeRange, queryRange, prob);
        String[] alignments = backtrack(biggestValue, dpTable, sequence, query, genomeRange, queryRange);
        this.genomeAlignment = alignments[0];
        this.queryAlignment = alignments[1];
        this.genomeRangeAignment = new int[]{genomeRange[0], genomeRange[0] + countGaps(genomeAlignment) - 1};
        this.queryRangeAlignment = new int[]{queryRange[0], queryRange[0] + countGaps(queryAlignment) - 1};

        StringBuilder sb = new StringBuilder();
        sb.append(sequence);
        dpTable = new ArrayList<>();
        setTableSize(genomeRange[1] - genomeRange[0] + 2, queryRange[1] - queryRange[0] + 2, dpTable);
        StringBuilder sb2 = new StringBuilder();
        sb2.append(query);
        int[] backwardsGenome = new int[2];
        int[] backwardsQuery = new int[2];
        backwardsGenome[1] = sequence.length() - 1 - genomeRange[0];
        backwardsGenome[0] = sequence.length() - 1 - genomeRange[1];
        backwardsQuery[1] = query.length() - 1 - queryRange[0];
        backwardsQuery[0] = query.length() - 1 - queryRange[1];
        List<Double> temp = Arrays.asList(prob);
        Collections.reverse(temp);
        Double[] reverseProb = (Double[])temp.toArray();

        biggestValue = SW(dpTable, sb.reverse().toString(), sb2.reverse().toString(), backwardsGenome, backwardsQuery, reverseProb);
        alignments = backtrack(biggestValue, dpTable, sb.toString(), sb2.toString(), backwardsGenome, backwardsQuery);
        if(alignments[0].length() > queryRange[1] - queryRange[0] + 1){
            StringBuilder newSeq = new StringBuilder();
            StringBuilder newQuer = new StringBuilder();
            int counter =  queryRange[1] - queryRange[0] + 1;
            while(counter < alignments[0].length()){
                newSeq.append(alignments[0].charAt(counter));
                newQuer.append(alignments[1].charAt(counter));
                counter++;
            }
            String newS = newSeq.toString();
            String newQ = newQuer.toString();
            this.genomeAlignment = newS + genomeAlignment;
            this.queryAlignment = newQ + queryAlignment;
            this.queryRangeAlignment[0] -= countGaps(newQ);
            this.genomeRangeAignment[0] -= countGaps(newS);

        }

        calculateProb(prob);

    }

    public int countGaps(String seq){
        int counter = 0;
        for(int i = 0; i < seq.length(); i++){
            if(seq.charAt(i) != '-'){
                counter++;
            }
        }
        return counter;
    }

    public void setTableSize(int seqDim, int queryDim, ArrayList<ArrayList<Pointer>> dpTable){
        while(dpTable.size() < seqDim){
            ArrayList<Pointer> column = new ArrayList<>();
            for(int i = 0; i < queryDim; i++){
                column.add(new Pointer());
            }
            dpTable.add(column);
        }
        for(int i = 0; i < seqDim; i++){
            for(int j = dpTable.get(i).size(); j < queryDim; j++){
                dpTable.get(i).add(new Pointer());
            }
        }
    }

    public Pointer SW(ArrayList<ArrayList<Pointer>> dpTable, String sequence, String query, int[] genomeRange, int[] queryRange, Double[] prob){

        //initialization
        setPosition(0,0,0.0, -1, -1, dpTable);
        for(int i = 0; i < dpTable.size(); i++){
            setPosition(i, 0, 0.0, -1, -1, dpTable);
        }
        for(int i = 0; i < dpTable.get(0).size(); i++){
            setPosition(0, i, 0.0, -1, -1, dpTable);
        }

        Pointer biggestValue = new Pointer();
        biggestValue.setPointer(0, 0, -1);

        int queryEnd = 0;
        int seqEnd = 0;
        while(queryEnd - biggestValue.prevIndexQuery < 10 && seqEnd - biggestValue.prevIndexSeq < 10){
            for(int sequenceIndex = 1; sequenceIndex < dpTable.size(); sequenceIndex++){
                for(int queryIndex = 1; queryIndex < dpTable.get(sequenceIndex).size(); queryIndex++){
                    if(queryIndex < queryEnd && sequenceIndex < seqEnd) {
                        continue;
                    }
                    double match = getValue(sequenceIndex-1, queryIndex-1, dpTable) +
                            getMatch(sequence, sequenceIndex - 1 + genomeRange[0], query.charAt(queryIndex - 1 + queryRange[0]), prob);
                    double deletion = getValue(sequenceIndex-1, queryIndex, dpTable) + getDeletion(sequenceIndex - 1 + genomeRange[0], prob);
                    double insertion = getValue(sequenceIndex, queryIndex-1, dpTable) + getDeletion(sequenceIndex - 1 + genomeRange[0], prob);
                    double[] maximum = max(match, deletion, insertion, 0);

                    int previousSeq = -1;
                    int previousQuery = -1;
                    if(maximum[0] == 0.0){
                        previousSeq = sequenceIndex - 1;
                        previousQuery = queryIndex - 1;
                    }else if(maximum[0] == 1.0){
                        previousSeq = sequenceIndex -1;
                        previousQuery = queryIndex;
                    }else if(maximum[0] == 2.0){
                        previousSeq = sequenceIndex;
                        previousQuery = queryIndex - 1;
                    }
                    setPosition(sequenceIndex, queryIndex, maximum[1], previousSeq, previousQuery, dpTable);
                    if(biggestValue.value < maximum[1]){
                        biggestValue.setPointer(sequenceIndex, queryIndex, maximum[1]);
                    }
                }
            }
            queryEnd = dpTable.size()-1;
            seqEnd = dpTable.size()-1;
            expandQuery(dpTable);
            expandSeq(dpTable);
            if(dpTable.size() - 1 + genomeRange[0] > sequence.length()){
                break;
            }
            if(dpTable.get(0).size() - 1 + queryRange[0] > query.length()){
                break;
            }
        }
        return biggestValue;
    }

    public String[] backtrack(Pointer biggestValue, ArrayList<ArrayList<Pointer>> dpTable, String sequence, String query, int[] genomeRange, int[] queryRange){
        int seqPos = biggestValue.prevIndexSeq;
        int querPos = biggestValue.prevIndexSeq;
        StringBuilder seqAlign = new StringBuilder();
        StringBuilder queryAlign = new StringBuilder();

        while(seqPos != -1 && querPos != -1){
            int nextSeq = dpTable.get(seqPos).get(querPos).prevIndexSeq;
            int nextQuer = dpTable.get(seqPos).get(querPos).prevIndexQuery;
            if(nextSeq == -1 || nextQuer == -1 || dpTable.get(nextSeq).get(nextQuer).value == 0){
                break;
            }
            int seqDif = seqPos - nextSeq;
            int querDif = querPos - nextQuer;

            if(seqDif == 1 && querDif == 1){
                seqAlign.append(sequence.charAt(seqPos -1 + genomeRange[0]));
                queryAlign.append(query.charAt(querPos - 1 + queryRange[0]));
            }else if(seqDif == 1){
                seqAlign.append(sequence.charAt(seqPos - 1 + genomeRange[0]));
                queryAlign.append("-");
            }else{
                seqAlign.append("-");
                queryAlign.append(query.charAt(querPos - 1 + queryRange[0]));
            }
            seqPos = nextSeq;
            querPos = nextQuer;

        }
        String[] value = {seqAlign.reverse().toString(), queryAlign.reverse().toString()};
        return value;
    }

    public void setPosition(int seqPos, int queryPos, Double value, int prevSeqPos, int prevQueryPos, ArrayList<ArrayList<Pointer>> dpTable){
        dpTable.get(seqPos).get(queryPos).setPointer(prevSeqPos, prevQueryPos, value);
    }

    public double getValue(int seqPos, int queryPos, ArrayList<ArrayList<Pointer>> dpTable){
        return dpTable.get(seqPos).get(queryPos).value;
    }

    public void expandSeq(ArrayList<ArrayList<Pointer>> dpTable){
        ArrayList<Pointer> column = new ArrayList<>();
        for(int i = 0; i < dpTable.get(0).size(); i++){
            column.add(new Pointer());
        }
        dpTable.add(column);
    }

    public void calculateProb(Double[] prob){
        int counter = 0;
        for(int i = 0; i < genomeAlignment.length(); i++){
            if(queryAlignment.charAt(i) == genomeAlignment.charAt(i)) {
                this.probability += prob[genomeRangeAignment[0] + counter];
                counter++;
            }else if(queryAlignment.charAt(i) != genomeAlignment.charAt(i) && genomeAlignment.charAt(i) != '-'){
                this.probability += prob[genomeRangeAignment[0] + counter];
                counter++;
            }
        }
        this.probability /= counter;
    }

    public String toString(){
        return (genomeAlignment + "\t" + genomeRangeAignment[0] + "-" + genomeRangeAignment[1] + "\n" + queryAlignment + "\t" +  queryRangeAlignment[0] + "-" + queryRangeAlignment[1] + "\n" + this.probability);
    }

    public void expandQuery(ArrayList<ArrayList<Pointer>> dpTable){
        for(int i = 0; i < dpTable.size(); i++){
           dpTable.get(i).add(new Pointer());
        }
    }

    public double getMatch(String sequence, int seqIndex, char query, Double[] probability){
        if(sequence.charAt(seqIndex) == query){
            return probability[seqIndex];
        }else{
            return -probability[seqIndex];
        }
    }

    //calculates the deletion penalty from index to index + 1
    public double getDeletion(int seqIndex, Double[] probability){
        return -probability[seqIndex];
    }

    public double[] max(double num1, double num2, double num3, double num4){
        double[] maximum = {-1, -Double.MAX_VALUE};
        if(num1 > maximum[1]){
            maximum = new double[]{0, num1};
        }
        if(num2 > maximum[1]){
            maximum = new double[]{1, num2};
        }
        if(num3 > maximum[1]){
            maximum = new double[]{2, num3};
        }
        if(num4 > maximum[1]){
            maximum = new double[]{3, num4};
        }
        return maximum;
    }

}

class Pointer{
    int prevIndexSeq;
    int prevIndexQuery;
    double value;

    public void setPointer(int prevIndexSeq, int prevIndexQuery, double value){
        this.prevIndexSeq = prevIndexSeq;
        this.prevIndexQuery = prevIndexQuery;
        this.value = value;
    }



}
