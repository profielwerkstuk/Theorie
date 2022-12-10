class Network {
	learningRate: number;
	limit: number;
	input: Layer;
	hidden: Layer[];
	output: Layer;
	constructor() {
		this.learningRate = 0.1;
		this.limit = 10000;
		this.input = new Layer(2, 0);
		this.hidden = [new Layer(3, 2), new Layer(3, 2), new Layer(3, 2)];
		this.output = new Layer(1, 3);
	}

	sigmoid(x: number) {
		return 1 / (1 + Math.exp(-x))
	}

	dsigmoid(x: number) {
		return x * (1 - x);
	}

	backPropagation(expectedValue: number) {
		for (let a = 0; a < this.output.nodes.length; a++) this.output.nodes[a].error = expectedValue - this.output.nodes[a].output;

		for (let a = 0; a < this.hidden.length; a++) {
			for (let b = 0; b < this.hidden[a].nodes.length; b++) {
				let total = 0;

				for (let c = 0; c < this.output.nodes.length; c++) total += this.output.nodes[c].error * this.output.nodes[c].weight[b];

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

	train(dataset: { input: any; output: any; }) {
		for (let n = 0; n < this.limit; n++) {
			for (let a = 0; a < dataset.input.length; a++) {
				this.predict(dataset.input[a]);
				this.backPropagation(dataset.output[a]);
			}
		}
	}

	predict(inputs: number[]) {
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
	nodes: Neuron[];
	constructor(totalNodes: number, totalWeights: number) {
		this.nodes = [];
		for (let a = 0; a < totalNodes; a++) this.nodes.push(new Neuron(totalWeights));
	}
}

class Neuron {
	weight: number[];
	error: number;
	output: number;
	constructor(totalWeights: number) {
		this.weight = [];
		this.error = 0;
		this.output = 1;

		for (let a = 0; a < totalWeights; a++) this.weight.push(Math.random() * 2 - 1);
	}
}

const dataSets = {
	"and": {
		input: [[0, 0], [1, 0], [0, 1], [1, 1]],
		output: [0, 0, 0, 1]
	},
	"or": {
		input: [[0, 0], [1, 0], [0, 1], [1, 1]],
		output: [0, 1, 1, 1]
	},
	"xor": {
		input: [[0, 0], [1, 0], [0, 1], [1, 1]],
		output: [0, 1, 1, 0]
	},
	"nand": {
		input: [[0, 0], [1, 0], [0, 1], [1, 1]],
		output: [1, 1, 1, 0]
	},
	"nor": {
		input: [[0, 0], [1, 0], [0, 1], [1, 1]],
		output: [1, 0, 0, 0]
	},
	"xnor": {
		input: [[0, 0], [1, 0], [0, 1], [1, 1]],
		output: [1, 0, 0, 1]
	}
}

const VerbumNetwork = new Network();
const dataset = dataSets.and;

console.log("Before training")
const before = getAccuracy(dataset, VerbumNetwork)
console.table(before);

let epoch = 1;
while (epoch--) VerbumNetwork.train(dataset);

console.log(`After training`)
const accuracy = getAccuracy(dataset, VerbumNetwork)
console.table(accuracy);

function getAccuracy(dataset: { [key: string]: any }, network: Network) {
	const table = [];
	for (const data in dataset.input) {
		const input = dataset.input[data];
		const expected = dataset.output[data];

		let accuracy = 0;
		for (let i = 0; i < 100; i++) {
			if (Math.round(VerbumNetwork.predict(input)) === expected) accuracy++
		}
		table.push({
			input: input.join(" "),
			expected: expected,
			accuracy: accuracy,
		})
	}

	return table;
}
