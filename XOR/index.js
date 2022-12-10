"use strict";
// input 0 0, output 0
// input 0 1, output 1
// input 1 0, output 1
// input 1 1, output 0
Object.defineProperty(exports, "__esModule", { value: true });
class Network {
    constructor() {
        this.learningRate = 0.1;
        this.limit = 10000;
        this.input = new Layer(2, 0);
        this.hidden = [new Layer(3, 2), new Layer(3, 2)];
        this.output = new Layer(1, 3);
    }
    sigmoid(x) {
        return 1 / (1 + Math.exp(-x));
    }
    dsigmoid(x) {
        return x * (1 - x);
    }
    backPropagation(expectedValue) {
        for (let a = 0; a < this.output.nodes.length; a++)
            this.output.nodes[a].error = expectedValue - this.output.nodes[a].output;
        for (let a = 0; a < this.hidden.length; a++) {
            for (let b = 0; b < this.hidden[a].nodes.length; b++) {
                let total = 0;
                for (let c = 0; c < this.output.nodes.length; c++)
                    total += this.output.nodes[c].error * this.output.nodes[c].weight[b];
                this.hidden[a].nodes[b].error = total;
            }
            for (let c = 0; c < this.hidden[a].nodes.length; c++) {
                for (let d = 0; d < this.input.nodes.length; d++) {
                    this.hidden[a].nodes[c].weight[d] += this.learningRate * this.hidden[a].nodes[c].error * this.input.nodes[d].output * this.dsigmoid(this.hidden[a].nodes[c].output);
                }
            }
        }
        for (let a = 0; a < this.output.nodes.length; a++) {
            for (let b = 0; b < this.hidden.length; b++) {
                for (let c = 0; c < this.hidden[b].nodes.length; c++) {
                    this.output.nodes[a].weight[c] += this.learningRate * this.output.nodes[a].error * this.hidden[b].nodes[c].output * this.dsigmoid(this.output.nodes[a].output);
                }
            }
        }
    }
    train(dataset) {
        for (let n = 0; n < this.limit; n++) {
            for (let a = 0; a < dataset.input.length; a++) {
                this.predict(dataset.input[a]);
                this.backPropagation(dataset.output[a]);
            }
        }
    }
    predict(inputs) {
        for (let i = 0; i < inputs.length; i++) {
            this.input.nodes[i].output = inputs[i];
        }
        for (let a = 0; a < this.hidden.length; a++) {
            for (let b = 0; b < this.hidden[a].nodes.length; b++) {
                let total = 0;
                for (let c = 0; c < this.input.nodes.length; c++)
                    total += this.input.nodes[c].output * this.hidden[a].nodes[b].weight[c];
                this.hidden[a].nodes[b].output = this.sigmoid(total);
            }
            for (let b = 0; b < this.output.nodes.length; b++) {
                let total = 0;
                for (let c = 0; c < this.hidden[a].nodes.length; c++)
                    total += this.hidden[a].nodes[c].output * this.output.nodes[b].weight[c];
                this.output.nodes[b].output = this.sigmoid(total);
            }
        }
        return this.output.nodes[0].output;
    }
}
class Layer {
    constructor(totalNodes, totalWeights) {
        this.nodes = [];
        for (let a = 0; a < totalNodes; a++)
            this.nodes.push(new Neuron(totalWeights));
    }
}
class Neuron {
    constructor(totalWeights) {
        this.weight = [];
        this.error = 0;
        this.output = 1;
        for (let a = 0; a < totalWeights; a++)
            this.weight.push(Math.random() * 2 - 1);
    }
}
const VerbumNetwork = new Network();
const dataset = {
    input: [[1, 1], [1, 0], [0, 0], [0, 1]],
    output: [0, 1, 0, 1]
};
VerbumNetwork.train(dataset);
console.log(Math.round(VerbumNetwork.predict([1, 0])));
console.log(Math.round(VerbumNetwork.predict([1, 1])));
console.log(Math.round(VerbumNetwork.predict([0, 0])));
console.log(Math.round(VerbumNetwork.predict([0, 1])));