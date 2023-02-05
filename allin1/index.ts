interface NEATConfig {
	populationSize: number;
	structure: StructureConfig;
	fitnessThreshold: number;
	maxEpoch: number;
	mutationRate?: MutationRateConfig;
	distanceConstants?: DistanceConfig;
	fitnessFunction: FitnessFunction;
}

interface MutationRateConfig {
	addNodeMR: number;
	addConnectionMR: number;
	removeNodeMR: number;
	removeConnectionMR: number;
	changeWeightMR: number;
}

interface DistanceConfig {
	c1: number;
	c2: number;
	c3: number;
	compatibilityThreshold: number;
}

interface FitnessFunction {
	(input: Genome): number;
}

interface StructureConfig {
	in: number;
	hidden: number;
	out: number;
	activationFunction: ActivationFunction;
}

interface ConnectionStructure {
	fNode: Neuron
	sNode: Neuron
}

interface ActivationFunction {
	(input: number, alpha?: number): number;
}

enum NodeType {
	INPUT = "INPUT",
	HIDDEN = "HIDDEN",
	OUTPUT = "OUTPUT"
}

const SIGMOID: ActivationFunction = (input: number): number => {
	return 1 / (1 + Math.exp(-input));
}

const TANH: ActivationFunction = (input: number): number => {
	if (input === Infinity) return 1;
	if (input === -Infinity) return -1;
	return (Math.exp(input) - Math.exp(-input)) / (Math.exp(input) + Math.exp(-input));
}

const RELU: ActivationFunction = (input: number): number => {
	return input > 0 ? input : 0;
}

const Gaussian: ActivationFunction = (input: number): number => {
	return Math.exp(-(input**2));
}

const expm1 = (x: number): number => {
	return Math.exp(x) - 1;
}

const ELU: ActivationFunction = (input: number, alpha = 0): number => {
	return ((input > 0) ? input : (alpha*expm1(input)));
}

const SELU: ActivationFunction = (input: number): number => {
	return 1.0507 * ELU(input, 1.67326);
}

class Neuron {

	value: number;
	innovation: number;
	type: NodeType;
	id: string;
	replacedConnection: Connection;
	active: boolean = true;
	inputCount: number = 0;
	inputTimes: number = 0;

	constructor(innovation: number, type: NodeType, replacedConnection?: Connection, id?: string, value?: number) {
		this.value = value ? value : 0;
		this.innovation = innovation;
		this.type = type;
		this.id = id ? id : this.newID();
		this.replacedConnection = replacedConnection ?? {} as Connection;
	}

	setValue(value: number) {
		this.value = value;
		if (this.type !== NodeType.INPUT) this.inputTimes++;
	}

	getValue(): number {
		return this.value;
	}

	getID(): string {
		return this.id;
	}

	getType(): NodeType {
		return this.type;
	}

	applyActivation(func: ActivationFunction) {
		this.value = func(this.value);
	}

	setNodeActivation(activation: boolean) {
		this.active = activation;
	}

	getState(): boolean {
		return this.inputTimes === this.inputCount;
	}

	newID(): string {
		const S4 = () => {
			return (((1 + Math.random()) * 65536) | 0).toString(16).substring(1);
		};
		return (`${S4()+S4()}-${S4()+S4()}`);
	}

	static getNodesByType(type: NodeType, nodes: Neuron[]): Neuron[] {
		let result: Neuron[] = [];
		for (let i = 0; i < nodes.length; i++) {
			if (nodes[i].getType() === type) result.push(nodes[i]);
		}
		return result;
	}

	static nodeExists(innovation: number, nodeDB: Neuron[]): number | undefined {
		for (let i = 0; i < nodeDB.length; i++) {
			if (nodeDB[i].replacedConnection.innovation === innovation) return nodeDB[i].innovation;
		}
		return undefined;
	}
}

class Connection {
	weight: number;
	active: boolean;
	input: Neuron;
	output: Neuron;
	innovation: number;

	constructor(input: Neuron, output: Neuron, innovation: number, weight?: number) {
		this.weight = weight ? weight : (Math.random() * 2) - 1;
		this.active = true;
		this.input = input;
		this.output = output;
		this.innovation = innovation;
	}

	randomizeWeight() {
		this.weight = (Math.random() * 2) - 1;
	}

	feedForward() {
		if (this.active) {
			this.output.setValue(this.output.getValue() + (this.input.getValue() * this.weight));
		}
	}

	getInputNode(): Neuron{
		return this.input;
	}

	getOutputNode(): Neuron{
		return this.output;
	}

	activateConnection() {
		this.active = true;
	}

	deactivateConnection() {
		this.active = false;
	}

	static isRecurrent(connection: Connection, genome: Genome): boolean {
		let node = connection.getInputNode();
		let stack = [connection];
		while (stack.length !== 0) {
			let connection = stack.shift();
			if (connection?.getOutputNode().getID() === node.getID()) return true;
			stack.push(
				...genome.connections.filter(gene => gene.getInputNode().getID() === connection?.getOutputNode().getID())
			);
		}
		return false;
	}

	static connectionExists(data: ConnectionStructure, connectionDB: Connection[]): number | undefined {
		for (let i = 0; i < connectionDB.length; i++) {
			if (data.fNode.innovation === connectionDB[i].getInputNode().innovation && data.sNode.innovation === connectionDB[i].getOutputNode().innovation) return connectionDB[i].innovation;
		}
		return undefined;
	}

	static inputConnectionsOfNode(node: Neuron, connections: Connection[]): Connection[] {
		let result: Connection[] = [];
		connections.forEach(connection => {
			if (connection.getInputNode().getID() === node.getID()) result.push(connection);
		});
		return result;
	}

	static outputConnectionsOfNode(node: Neuron, connections: Connection[]): Connection[] {
		let result: Connection[] = [];
		connections.forEach(connection => {
			if (connection.getOutputNode().getID() === node.getID()) result.push(connection);
		});
		return result;
	}
}

class Genome {
	nodes: Neuron[] = [];
	connections: Connection[] = [];
	fitness: number = 0;
	config: StructureConfig;
	activationFunction: ActivationFunction;

	constructor(config: StructureConfig) {
		this.activationFunction = config.activationFunction;
		this.config = config;
		for (let i = 0; i < config.in; i++) {
			this.nodes.push(new Neuron(i, NodeType.INPUT));
		}

		for (let i = config.in; i < config.in + config.hidden; i++) {
			this.nodes.push(new Neuron(i, NodeType.HIDDEN));
		}

		for (let i = config.in + config.hidden; i < config.in + config.hidden + config.out; i++) {
			this.nodes.push(new Neuron(i, NodeType.OUTPUT));
		}
	}

	activate(input: number[]): number[] {
		for (let i = 0; i < this.nodes.length; i++) {
			this.nodes[i].inputCount = Connection.outputConnectionsOfNode(this.nodes[i], this.connections).length;
			this.nodes[i].inputTimes = 0;
			this.nodes[i].value = 0;
		}

		let stack =Neuron.getNodesByType(NodeType.INPUT, this.nodes);
		stack = stack.sort((a, b) => a.innovation < b.innovation ? -1 : 1);
		for (let i = 0; i < stack.length; i++) {
			stack[i].setValue(input[i]);
		}

		while (stack.length) {
			let node = stack.splice(stack.indexOf(stack.filter(n => n.getState())[0]), 1)[0];
			let connections = Connection.inputConnectionsOfNode(node, this.connections);
			connections.forEach(connection => {
				connection.feedForward();
				if (connection.getOutputNode().getState()) {
					connection.getOutputNode().inputTimes = 0;
					connection.getOutputNode().applyActivation(this.activationFunction);
					stack.push(connection.getOutputNode());
				}
			});
		}

		return this.getOutputValues();
	}

	getNodes() {
		return this.nodes;
	}

	getOutputValues(): number[] {
		let oNodes =Neuron.getNodesByType(NodeType.OUTPUT, this.nodes);
		oNodes = oNodes.sort((a, b) => a.innovation < b.innovation ? -1 : 1);
		let result: number[] = [];
		oNodes.forEach(node => {
			result.push(node.getValue());
		});
		return result;
	}

	hasConnection(innovation: number): Connection | boolean {
		for (let i = 0; i < this.connections.length; i++) {
			if (this.connections[i].innovation === innovation) return this.connections[i];
		}
		return false;
	}

	hasNode(innovation: number): Neuron| boolean {
		for (let i = 0; i < this.nodes.length; i++) {
			if (this.nodes[i].innovation === innovation) return this.nodes[i];
		}
		return false;
	}

	randomConnectionStructure(): ConnectionStructure | void {
		let tries = 0;
		let fNode = this.nodes[Math.floor(Math.random() * this.nodes.length)];
		let sNode = this.nodes[Math.floor(Math.random() * this.nodes.length)];
		while (fNode.id === sNode.id || (fNode.getType() === NodeType.INPUT && sNode.getType() === NodeType.INPUT) || (fNode.getType() === NodeType.OUTPUT && sNode.getType() === NodeType.OUTPUT)) {
			sNode = this.nodes[Math.floor(Math.random() * this.nodes.length)];
			tries++;
		}
		if (!(tries > 20 || fNode.getType() === NodeType.OUTPUT || sNode.getType() === NodeType.INPUT)) return { fNode: fNode, sNode: sNode };
		else return;
	}

	addNode(rConnection: Connection, neat: NEAT): Neuron{
		let nInnovation =Neuron.nodeExists(rConnection.innovation, neat.nodeDB);

		if (nInnovation) {
			let existing = this.hasNode(nInnovation);
			if (!existing) {
				let newNode = new Neuron(nInnovation, NodeType.HIDDEN, rConnection);
				this.nodes.push(newNode);
				return newNode;
			} else {
				// @ts-ignore
				existing.setNodeActivation(true);
				// @ts-ignore
				return existing;
			}
		} else {
			neat.nodeInnovation++;
			let newNode = new Neuron(neat.nodeInnovation, NodeType.HIDDEN, rConnection);
			this.nodes.push(newNode);
			neat.nodeDB.push(newNode);
			// @ts-ignore
			return newNode;
		}
	}

	addConnection(tNodes: ConnectionStructure, neat: NEAT): Connection | void {
		let innovation = Connection.connectionExists(tNodes, neat.connectionDB);
		if (innovation) {
			let existing = this.hasConnection(innovation);
			if (!existing) {
				let newConnection = new Connection(tNodes.fNode, tNodes.sNode, innovation);
				if (Connection.isRecurrent(newConnection, this)) {
					return;
				} else {
					this.connections.push(newConnection);
					return newConnection;
				}
			}
		} else {
			neat.connectionInnovation++;
			let newConnection = new Connection(tNodes.fNode, tNodes.sNode, neat.connectionInnovation);
			if (!Connection.isRecurrent(newConnection, this)) {
				neat.connectionDB.push(newConnection);
				this.connections.push(newConnection);
				return newConnection;
			} else {
				neat.connectionInnovation--;
				return;
			}
		}
	}

	getGenesByInnovation(): Connection[] {
		let result: Connection[] = [];
		for (let i = 0; i < this.connections.length; i++) {
			result[this.connections[i].innovation] = this.connections[i];
		}
		return result;
	}

	mutateWeights(rate: number) {
		for (let i = 0; i < this.connections.length; i++) {
			if (Math.random() < rate) {
				this.connections[i].randomizeWeight();
			}
		}
	}

	mutateConnection(neat: NEAT) {
		let tNodes = this.randomConnectionStructure();
		if (tNodes) this.addConnection(tNodes, neat);
	}

	mutateDeactivateConnection() {
		let rConnection = this.connections[Math.floor(Math.random() * this.connections.length)];
		if (rConnection) rConnection.deactivateConnection();
	}

	mutateNode(neat: NEAT) {
		let rConnection = this.connections[Math.floor(Math.random() * this.connections.length)];

		if (rConnection) {
			if (!rConnection.active) return;
			rConnection.deactivateConnection();
			let iNode = rConnection.getInputNode();
			let oNode = rConnection.getOutputNode();

			let node = this.addNode(rConnection, neat);
			let fConnection = { fNode: iNode, sNode: node };
			let sConnection = { fNode: node, sNode: oNode };
			this.addConnection(fConnection, neat);
			this.addConnection(sConnection, neat);
		}
	}

	mutateDeactivateNode() {
		let node = this.nodes[Math.floor(Math.random() * this.nodes.length)];
		if (node.replacedConnection) {
			node.setNodeActivation(false);
			for (let i = 0; i < this.connections.length; i++) {
				if (this.connections[i].getInputNode().getID() === node.getID() || this.connections[i].getOutputNode().getID() === node.getID()) this.connections[i].deactivateConnection();
			}
		}
	}

	addGene(gene: Connection) {
		let iNode = gene.getInputNode();
		let oNode = gene.getOutputNode();

		let childiNode = this.hasNode(iNode.innovation);
		let childoNode = this.hasNode(oNode.innovation);

		let iNodeConnection;
		let oNodeConnection;
		if (!childiNode) {
			iNodeConnection = new Neuron(iNode.innovation, iNode.getType(), iNode.replacedConnection);
			this.nodes.push(iNodeConnection);
		} else {
			iNodeConnection = childiNode;
			// @ts-ignore
			childiNode.setNodeActivation(true);
		}

		if (!childoNode) {
			oNodeConnection = new Neuron(oNode.innovation, oNode.getType(), oNode.replacedConnection);
			this.nodes.push(oNodeConnection);
		} else {
			oNodeConnection = childoNode;
			// @ts-ignore
			childoNode.setNodeActivation(true);
		}

		let childConnection = this.hasConnection(gene.innovation);
		if (!childConnection && iNodeConnection != true && oNodeConnection != true) {
			let connection = new Connection(iNodeConnection, oNodeConnection, gene.innovation, gene.weight);
			if (!Connection.isRecurrent(connection, this)) this.connections.push(new Connection(iNodeConnection, oNodeConnection, gene.innovation, gene.weight));
		} else {
			// @ts-ignore
			childConnection.activateConnection();
		}
	}

	static crossover(genome1: Genome, genome2: Genome, config: StructureConfig): Genome {
		let child = new Genome(config);
		const [hFit, lFit] = [genome1, genome2].sort((a, b) => b.fitness - a.fitness);
		const hFitGenes = hFit.getGenesByInnovation();
		const lFitGenes = lFit.getGenesByInnovation();

		for (let i = 0; i < Math.max(hFitGenes.length, lFitGenes.length); i++) {
			if (hFitGenes[i] !== undefined && lFitGenes[i] !== undefined) {
				if (Math.random() < 0.5) {
					child.addGene(hFitGenes[i]);
				} else {
					child.addGene(lFitGenes[i]);
				}
			} else if (hFitGenes[i] !== undefined) {
				child.addGene(hFitGenes[i]);
			} else if (lFitGenes[i] !== undefined) {
				child.addGene(lFitGenes[i]);
			}
		}
		return child;
	}

	// Using variable names used in the original paper.
	static isCompatible(genome1: Genome, genome2: Genome, config: DistanceConfig): boolean {
		let genes1 = genome1.getGenesByInnovation();
		let genes2 = genome2.getGenesByInnovation();
		let E = Math.abs(genes1.length - genes2.length);
		let N = (Math.max(genes1.length, genes2.length) < 20) ? 1 : Math.max(genes1.length, genes2.length);
		let D = 0;

		for (let i = 0; i < Math.min(genes1.length, genes2.length); i++) {
			if ((genes1[i] === undefined && genes2[i] !== undefined) || (genes1[i] !== undefined && genes2[i] === undefined)) D++;
		}

		let W = 0;
		let count = 0;
		for (let i = 0; i < Math.min(genes1.length, genes2.length); i++) {
			if (genes1[i] !== undefined && genes2[i] !== undefined) {
				W += Math.abs(genes1[i].weight - genes2[i].weight);
				count++;
			}
		}
		W /= count;
		W = isNaN(W) ? 0 : W;
		return (((config.c1 * E) / N) + ((config.c2 * D) / N) + config.c3 * W) < config.compatibilityThreshold;
	}
}

class Species {
	genomes: Genome[] = [];
	populationCap: number = 0;
	adjustedFitness: number = 0;

	addGenome(genome: Genome) {
		this.genomes.push(genome);
	}

	getGenomes(): Genome[] {
		return this.genomes;
	}

	getSpecimen(): Genome {
		return this.genomes[0];
	}

	randomGenome(): Genome {
		return this.genomes[Math.floor(Math.random() * this.genomes.length)];
	}

	adjustFitness(): number {
		this.adjustedFitness = 0;
		this.genomes.forEach(genome => {
			this.adjustedFitness += genome.fitness / this.genomes.length;
		});
		return this.adjustedFitness;
	}

	repopulate(config: StructureConfig) {
		this.genomes = this.genomes.sort((a, b) => b.fitness - a.fitness)
		let half_length = Math.ceil(this.genomes.length / 2);
		this.genomes = this.genomes.splice(0, half_length);

		let newGenomes: Genome[] = [];
		while (newGenomes.length < this.populationCap) {
			newGenomes.push(Genome.crossover(this.randomGenome(), this.randomGenome(), config));
		}
		this.genomes = newGenomes;
	}

	mutateConnection(rate: number, neat: NEAT) {
		for (let i = 0; i < this.genomes.length; i++) {
			if (Math.random() < rate) {
				this.genomes[i].mutateConnection(neat);
			}
		}
	}

	mutateDeactivateConnection(rate: number) {
		for (let i = 0; i < this.genomes.length; i++) {
			if (Math.random() < rate) {
				this.genomes[i].mutateDeactivateConnection();
			}
		}
	}

	mutateNode(rate: number, neat: NEAT) {
		for (let i = 0; i < this.genomes.length; i++) {
			if (Math.random() < rate) {
				this.genomes[i].mutateNode(neat);
			}
		}
	}

	mutateDeactivateNode(rate: number) {
		for (let i = 0; i < this.genomes.length; i++) {
			if (Math.random() < rate) {
				this.genomes[i].mutateDeactivateNode();
			}
		}
	}

	mutateWeight(rate: number) {
		for (let i = 0; i < this.genomes.length; i++) {
			this.genomes[i].mutateWeights(rate);
		}
	}

	static speciate(genomes: Genome[], config: DistanceConfig): Species[] {
		let species: Species[] = [];
		let firstSpecies = new Species();
		firstSpecies.addGenome(genomes[0]);
		species.push(firstSpecies);
		for (let i = 1; i < genomes.length; i++) {
			let foundMatch = false;
			for (let o = 0; o < species.length; o++) {
				if (Genome.isCompatible(genomes[i], species[o].getSpecimen(), config)) {
					species[o].addGenome(genomes[i]);
					foundMatch = true;
					break;
				}
			}
			if (!foundMatch) {
				let newSpecies = new Species();
				newSpecies.addGenome(genomes[i]);
				species.push(newSpecies);
			}
		}
		return species;
	}
}

class NEAT {
	config: NEATConfig;
	species: Species[] = [];
	nodeInnovation: number;
	connectionInnovation: number;
	connectionDB: Connection[] = [];
	nodeDB: Neuron[] = [];
	epoch: number = 0;

	constructor(config: NEATConfig) {

		this.config = config;

		config.structure.hidden = (config.structure.hidden !== undefined) ? config.structure.hidden : 0;
		this.nodeInnovation = config.structure.in + config.structure.hidden + config.structure.out - 1;
		this.connectionInnovation = 0;

		if (config.mutationRate) {
			config.mutationRate.addNodeMR = (config.mutationRate.addNodeMR !== undefined) ? config.mutationRate.addNodeMR : 0.01;
			config.mutationRate.addConnectionMR = (config.mutationRate.addConnectionMR !== undefined) ? config.mutationRate.addConnectionMR : 0.02;
			config.mutationRate.removeNodeMR = (config.mutationRate.removeNodeMR !== undefined) ? config.mutationRate.removeNodeMR : 0.005;
			config.mutationRate.removeConnectionMR = (config.mutationRate.removeConnectionMR !== undefined) ? config.mutationRate.removeConnectionMR : 0.005;
			config.mutationRate.changeWeightMR = (config.mutationRate.changeWeightMR !== undefined) ? config.mutationRate.changeWeightMR : 0.01;
		} else {
			config.mutationRate = { addNodeMR: 0.01, addConnectionMR: 0.02, removeNodeMR: 0.005, removeConnectionMR: 0.005, changeWeightMR: 0.01 };
		}

		this.species.push(new Species());
		for (let i = 0; i < config.populationSize; i++) {
			this.species[0].addGenome(new Genome({ in: config.structure.in, hidden: config.structure.hidden, out: config.structure.out, activationFunction: config.structure.activationFunction }));
		}
	}

	mutate() {
		this.species.forEach(specie => {
			specie.mutateNode(this.config.mutationRate!.addNodeMR, this);
			specie.mutateConnection(this.config.mutationRate!.addConnectionMR, this);
			specie.mutateDeactivateNode(this.config.mutationRate!.removeNodeMR);
			specie.mutateDeactivateConnection(this.config.mutationRate!.removeConnectionMR);
			specie.mutateWeight(this.config.mutationRate!.changeWeightMR);
		});
	}

	speciate() {
		let genomes: Genome[] = [];
		for (let i = 0; i < this.species.length; i++) {
			genomes = genomes.concat(this.species[i].getGenomes());
		}

		this.species = Species.speciate(genomes, this.config.distanceConstants!);
	}

	assignPopulationLimit() {
		let total = 0;
		this.species.forEach(specie => {
			total += specie.adjustFitness();
		});

		this.species.forEach(specie => {
			let normalized = specie.adjustedFitness / total;
			specie.populationCap = Math.floor(normalized * this.config.populationSize);
			if (isNaN(specie.populationCap) || specie.populationCap < 0) {
				specie.populationCap = 0;
			}
		});

		for (let i = this.species.length - 1; i >= 0; i--) {
			if (this.species[i].populationCap === 0 || this.species[i].genomes.length < 1) this.species.splice(i, 1);
		}
	}

	repopulate() {
		this.species.forEach(specie => {
			specie.repopulate(this.config.structure);
		});
	}

	run(): Genome[] | undefined {
		let fitness: Genome[];
		while (this.config.maxEpoch > this.epoch) {
			fitness = [];
			let genomes: Genome[] = [];
			for (let i = 0; i < this.species.length; i++) {
				genomes = genomes.concat(this.species[i].getGenomes());
			}

			genomes.forEach(genome => {
				genome.fitness = this.config.fitnessFunction(genome);
				fitness.push(genome);
				if (isNaN(genome.fitness) || genome.fitness === undefined) throw new Error("Fitness function returned NaN or undefined.");
			});

			if (fitness.filter(genome => genome.fitness > this.config.fitnessThreshold).length > 0) return fitness.filter(genome => genome.fitness > this.config.fitnessThreshold);

			this.speciate();
			this.assignPopulationLimit();
			this.repopulate();
			this.mutate();
			this.epoch++;
		}
	}
}

let best: {score: number, genome: any} = {
    score: 0,
    genome: null,
};

function fitnessFunction(a: any) {
    let fitness = 4;
    fitness -= Math.abs(a.activate([1, 1])[0]);
    fitness -= Math.abs(a.activate([1, 0])[0] - 1);
    fitness -= Math.abs(a.activate([0, 1])[0] - 1);
    fitness -= Math.abs(a.activate([0, 0])[0]);
    if (a.connections.length < 2) fitness *= 0.001;

    const score = Math.max(fitness, 0.001)

    if (score > best.score) {
        best.score = score;
        best.genome = a;
    }

    return score;
};

let config = {
    populationSize: 9999,
    structure: {
        in: 2,
        hidden: 0,
        out: 1,
        activationFunction: RELU
    },
    mutationRate: {
        addNodeMR: 0.005,
        addConnectionMR: 0.01,
        removeNodeMR: 0.0001,
        removeConnectionMR: 0.01,
        changeWeightMR: 0.1
    },
    distanceConstants: {
        c1: 2,
        c2: 0.5,
        c3: 1,
        compatibilityThreshold: 1.5
    },
    fitnessThreshold: 3.5,
    fitnessFunction: fitnessFunction,
    maxEpoch: Math.round(Math.exp(50)),
};

const neat = new NEAT(config);
console.log('Starting...');

neat.run();

console.log('Best score: ' + best.score);
console.log('Best genome: ' + JSON.stringify(best.genome));
console.log("Testing best genome");
console.log(0, ":", Math.round(best.genome.activate([1, 1])[0]), Math.round(best.genome.activate([0, 0])[0]));
console.log(1, ":", Math.round(best.genome.activate([1, 0])[0]), Math.round(best.genome.activate([0, 1])[0]));

